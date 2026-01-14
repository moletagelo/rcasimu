function [rf_RC, tstart_RC, rf_CR, tstart_CR] = acquire_rf_data(...
    Tx_row, Rx_col, Tx_col, Rx_row, positions, amplitudes, params)
%ACQUIRE_RF_DATA 采集Row-Column RF数据 - 稳定并行版
%   修复了频繁初始化导致的 Field II 崩溃问题

    N = params.N_elements;
    
    %% 准备并行环境
    fprintf('    正在配置并行计算资源...\n');
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        % 如果没有并行池，启动一个（根据电脑核数自动决定数量）
        parpool('local'); 
    end
    
    % 预先定义好所有 worker 需要的参数，避免传递大对象
    fs = params.fs;
    c = params.c;
    f0 = params.f0;
    lambda = params.lambda;
    element_size = params.element_size;
    kerf = params.kerf;

    %% Row发射 -> Column接收
    rf_RC = cell(N, 1);
    tstart_RC = zeros(N, 1);
    
    fprintf('    开始并行计算 (Row Tx)...\n');
    
    parfor tx = 1:N
        % 1. 确保当前 Worker 已初始化 Field II (关键修复)
        ensure_field_init(fs, c);
        
        % 2. 在 Worker 本地重新创建换能器
        % (Field II 的对象 ID 是本地的，不能跨 Worker 共享)
        [Tx_local, Rx_local] = create_local_row_col(N, element_size, kerf, f0, tx, 1);
        
        % 3. 执行计算 (带错误保护)
        try
            [rf, ts] = calc_scat_batched(Tx_local, Rx_local, positions, amplitudes, fs);
            if ~isempty(rf) && numel(rf) > 50
                rf_RC{tx} = single(rf);
                tstart_RC(tx) = ts;
            end
        catch ME
            fprintf('    Row Tx %d 出错: %s\n', tx, ME.message);
        end
        
        % 4. 释放本地换能器内存 (但不关闭 Field II)
        xdc_free(Tx_local);
        xdc_free(Rx_local);
    end
    
    %% Column发射 -> Row接收
    rf_CR = cell(N, 1);
    tstart_CR = zeros(N, 1);
    
    fprintf('    开始并行计算 (Col Tx)...\n');
    
    parfor tx = 1:N
        ensure_field_init(fs, c);
        
        % mode=2 表示 Col 发射, Row 接收
        [Tx_local, Rx_local] = create_local_row_col(N, element_size, kerf, f0, tx, 2);
        
        try
            [rf, ts] = calc_scat_batched(Tx_local, Rx_local, positions, amplitudes, fs);
            if ~isempty(rf) && numel(rf) > 50
                rf_CR{tx} = single(rf);
                tstart_CR(tx) = ts;
            end
        catch ME
            fprintf('    Col Tx %d 出错: %s\n', tx, ME.message);
        end
        
        xdc_free(Tx_local);
        xdc_free(Rx_local);
    end
    
    fprintf('    有效数据: RC=%d, CR=%d\n', ...
            sum(~cellfun(@isempty, rf_RC)), sum(~cellfun(@isempty, rf_CR)));
end

%% === 辅助函数 ===

function ensure_field_init(fs, c)
%ENSURE_FIELD_INIT 确保当前 Worker 已初始化 Field II
%   使用 persistent 变量避免重复初始化
    persistent is_initialized
    if isempty(is_initialized) || ~is_initialized
        % 尝试一个简单命令测试是否存活
        try
            set_field('c', c);
        catch
            % 如果失败或未初始化，则重新初始化
            field_init(-1);
            set_field('fs', fs);
            set_field('c', c);
            set_field('use_triangles', 0); % 关闭三角形加速以提高稳定性
        end
        is_initialized = true;
    end
end

function [Tx, Rx] = create_local_row_col(N, width, kerf, f0, tx_idx, mode)
%CREATE_LOCAL_ROW_COL 在 Worker 内部创建换能器配置
    
    height = 5e-3; % 假设高度
    
    % 创建线性阵列
    % 注意：Field II 对象是指针，必须在同一个线程内创建和使用
    Tx_array = xdc_linear_array(N, width, height, kerf, 1, 1, [0 0 Inf]);
    Rx_array = xdc_linear_array(N, width, height, kerf, 1, 1, [0 0 Inf]);
    
    % 设置脉冲响应
    impulse_resp = gauspuls('cutoff', f0, 0.6, [], -6);
    xdc_impulse(Tx_array, impulse_resp);
    xdc_impulse(Rx_array, impulse_resp);
    xdc_excitation(Tx_array, impulse_resp);
    
    % 设置聚焦和变迹
    if mode == 1 % Row -> Col
        Tx = Tx_array;
        Rx = Rx_array;
    else % Col -> Row
        Tx = Rx_array; % 交换物理意义
        Rx = Tx_array;
    end
    
    % 发射：单阵元激活
    apo_tx = zeros(1, N);
    apo_tx(tx_idx) = 1;
    xdc_apodization(Tx, 0, apo_tx);
    xdc_focus(Tx, 0, [0 0 100]); % 远场聚焦（平面波）
    
    % 接收：全孔径
    xdc_apodization(Rx, 0, ones(1, N));
    xdc_focus(Rx, 0, [0 0 0]); % 动态聚焦通常在波束形成中做，这里设0
end

function [rf_total, tstart_total] = calc_scat_batched(Tx, Rx, pos, amp, fs)
%CALC_SCAT_BATCHED 分批计算散射
    BATCH_SIZE = 500; 
    n_points = size(pos, 1);
    
    % 1. 小数据量直接计算 (加了 try-catch 保护)
    if n_points <= BATCH_SIZE
        try
            [rf_total, tstart_total] = calc_scat(Tx, Rx, pos, amp);
        catch
            rf_total = [];
            tstart_total = 0;
        end
        return;
    end
    
    % 2. 大数据量分批
    n_batches = ceil(n_points / BATCH_SIZE);
    rf_parts = cell(n_batches, 1);
    tstart_parts = zeros(n_batches, 1);
    valid_mask = false(n_batches, 1);
    
    for b = 1:n_batches
        idx_start = (b-1)*BATCH_SIZE + 1;
        idx_end = min(b*BATCH_SIZE, n_points);
        
        try
            [rf, ts] = calc_scat(Tx, Rx, pos(idx_start:idx_end, :), amp(idx_start:idx_end));
            if ~isempty(rf)
                rf_parts{b} = rf;
                tstart_parts(b) = ts;
                valid_mask(b) = true;
            end
        catch
            % 单个批次失败不影响整体
        end
    end
    
    if ~any(valid_mask)
        rf_total = [];
        tstart_total = 0;
        return;
    end
    
    % 合并逻辑
    valid_indices = find(valid_mask);
    min_tstart = min(tstart_parts(valid_indices));
    tstart_total = min_tstart;
    
    % 计算最大样本长度
    max_len = 0;
    for i = 1:length(valid_indices)
        idx = valid_indices(i);
        offset = round((tstart_parts(idx) - min_tstart) * fs);
        end_idx = offset + size(rf_parts{idx}, 1);
        if end_idx > max_len, max_len = end_idx; end
    end
    
    rf_total = zeros(max_len, size(rf_parts{valid_indices(1)}, 2));
    
    for i = 1:length(valid_indices)
        idx = valid_indices(i);
        offset = round((tstart_parts(idx) - min_tstart) * fs);
        current_rf = rf_parts{idx};
        s = offset + 1;
        e = offset + size(current_rf, 1);
        rf_total(s:e, :) = rf_total(s:e, :) + current_rf;
    end
end
