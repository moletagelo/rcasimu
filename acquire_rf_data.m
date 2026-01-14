function [rf_RC, tstart_RC, rf_CR, tstart_CR] = acquire_rf_data(...
    Tx_row, Rx_col, Tx_col, Rx_row, positions, amplitudes, params)
%ACQUIRE_RF_DATA 采集Row-Column RF数据 - 修复几何定义的并行版

    N = params.N_elements;
    
    %% 准备并行环境
    fprintf('    正在配置并行计算资源...\n');
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local'); 
    end
    
    % 提取所有需要的参数，避免传递 params 结构体
    fs = params.fs;
    c = params.c;
    f0 = params.f0;
    width = params.element_size;
    height_full = params.array_size; % 阵元全长
    kerf = params.kerf;
    row_pos_x = params.row_pos_x;
    col_pos_y = params.col_pos_y;
    
    %% Row发射 -> Column接收
    rf_RC = cell(N, 1);
    tstart_RC = zeros(N, 1);
    
    fprintf('    开始并行计算 (Row Tx)...\n');
    
    parfor tx = 1:N
        ensure_field_init(fs, c);
        
        % mode=1: Row发射(沿X排列, 长边沿Y), Col接收(沿Y排列, 长边沿X)
        [Tx_local, Rx_local] = create_local_transducers(N, width, height_full, kerf, f0, row_pos_x, col_pos_y, tx, 1);
        
        try
            [rf, ts] = calc_scat_batched(Tx_local, Rx_local, positions, amplitudes, fs);
            if ~isempty(rf) && numel(rf) > 50
                rf_RC{tx} = single(rf);
                tstart_RC(tx) = ts;
            end
        catch ME
            fprintf('    Row Tx %d 出错: %s\n', tx, ME.message);
        end
        
        xdc_free(Tx_local);
        xdc_free(Rx_local);
    end
    
    %% Column发射 -> Row接收
    rf_CR = cell(N, 1);
    tstart_CR = zeros(N, 1);
    
    fprintf('    开始并行计算 (Col Tx)...\n');
    
    parfor tx = 1:N
        ensure_field_init(fs, c);
        
        % mode=2: Col发射, Row接收
        [Tx_local, Rx_local] = create_local_transducers(N, width, height_full, kerf, f0, row_pos_x, col_pos_y, tx, 2);
        
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
    persistent is_initialized
    if isempty(is_initialized) || ~is_initialized
        try
            set_field('c', c);
        catch
            field_init(-1);
            set_field('fs', fs);
            set_field('c', c);
            set_field('use_triangles', 0);
        end
        is_initialized = true;
    end
end

function [Tx, Rx] = create_local_transducers(N, el_width, array_len, kerf, f0, row_pos_x, col_pos_y, tx_idx, mode)
%CREATE_LOCAL_TRANSDUCERS 在Worker中重建正确的几何
    
    % --- 脉冲响应 ---
    impulse_resp = sin(2*pi*f0*(-2/f0:1/f0/10:2/f0));
    impulse_resp = impulse_resp .* hanning(length(impulse_resp))';
    
    % --- 定义 Row 阵列 (沿X排列，阵元长边沿Y) ---
    % 使用 xdc_linear_array，它默认创建沿X轴排列的阵列
    % 这里的 width 是 x方向尺寸，height 是 y方向尺寸(即array_len)
    Row_Array = xdc_linear_array(N, el_width, array_len, kerf, 1, 4, [0 0 50e-3]); 
    xdc_impulse(Row_Array, impulse_resp);
    
    % --- 定义 Column 阵列 (沿Y排列，阵元长边沿X) ---
    % 必须手动构建矩形，因为 xdc_linear_array 只能沿X轴排列
    rect_col = zeros(N, 19);
    for i = 1:N
        rect_col(i, 1) = i;
        rect_col(i, 2) = array_len;      % X尺寸 (长边)
        rect_col(i, 3) = el_width;       % Y尺寸 (短边)
        rect_col(i, 4:6) = [4 1 1];      % 子阵元划分 (X方向分4份以提高精度)
        rect_col(i, 7) = 0;              % X中心
        rect_col(i, 8) = col_pos_y(i);   % Y中心 (根据params分布)
        rect_col(i, 9) = 0;              % Z中心
        rect_col(i, 11:19) = [0 0 1 1 0 0 0 0 0]; % 欧拉角和标志位
    end
    Col_Array = xdc_rectangles(rect_col, [0 0 50e-3], [0 0 1e9]);
    xdc_impulse(Col_Array, impulse_resp);
    
    % --- 配置发射/接收 ---
    if mode == 1 
        % Case: Row发射 -> Col接收
        Tx = Row_Array;
        Rx = Col_Array;
        
        % 发射聚焦
        apo = zeros(1, N);
        apo(tx_idx) = 1;
        xdc_apodization(Tx, 0, apo);
        xdc_excitation(Tx, impulse_resp); % 简单起见用脉冲响应作为激励
        xdc_focus(Tx, 0, [0 0 100]);      % 平面波
        
        % 接收设置
        xdc_apodization(Rx, 0, ones(1, N));
        xdc_focus(Rx, 0, [0 0 0]);
        
    else 
        % Case: Col发射 -> Row接收
        Tx = Col_Array;
        Rx = Row_Array;
        
        apo = zeros(1, N);
        apo(tx_idx) = 1;
        xdc_apodization(Tx, 0, apo);
        xdc_excitation(Tx, impulse_resp);
        xdc_focus(Tx, 0, [0 0 100]);
        
        xdc_apodization(Rx, 0, ones(1, N));
        xdc_focus(Rx, 0, [0 0 0]);
    end
end

function [rf_total, tstart_total] = calc_scat_batched(Tx, Rx, pos, amp, fs)
    BATCH_SIZE = 500; 
    n_points = size(pos, 1);
    
    if n_points <= BATCH_SIZE
        try
            [rf_total, tstart_total] = calc_scat(Tx, Rx, pos, amp);
        catch
            rf_total = []; tstart_total = 0;
        end
        return;
    end
    
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
        end
    end
    
    if ~any(valid_mask)
        rf_total = []; tstart_total = 0; return;
    end
    
    valid_indices = find(valid_mask);
    min_tstart = min(tstart_parts(valid_indices));
    tstart_total = min_tstart;
    
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
        rf = rf_parts{idx};
        rf_total(offset+1 : offset+size(rf,1), :) = rf_total(offset+1 : offset+size(rf,1), :) + rf;
    end
end
