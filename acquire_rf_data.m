function [rf_RC, tstart_RC, rf_CR, tstart_CR] = acquire_rf_data(...
    Tx_row, Rx_col, Tx_col, Rx_row, positions, amplitudes, params)
%ACQUIRE_RF_DATA 采集Row-Column RF数据
%   分批次采集以节省内存 - 修复了Field II内存不足的问题

    N = params.N_elements;
    
    %% Row发射 -> Column接收
    rf_RC = cell(N, 1);
    tstart_RC = zeros(N, 1);
    
    for tx = 1:N
        % 设置发射孔径(单阵元)
        apo_tx = zeros(1, N);
        apo_tx(tx) = 1;
        xdc_apodization(Tx_row, 0, apo_tx);
        xdc_focus(Tx_row, 0, [0 0 500e-3]);  % 平面波
        
        % 设置接收孔径(全孔径)
        xdc_apodization(Rx_col, 0, ones(1, N));
        
        % 计算散射响应 (使用分批处理函数)
        try
            [rf, ts] = calc_scat_batched(Tx_row, Rx_col, positions, amplitudes, params.fs);
            
            if ~isempty(rf) && numel(rf) > 50
                rf_RC{tx} = single(rf);  % 使用单精度节省内存
                tstart_RC(tx) = ts;
            end
        catch ME
            fprintf('    警告: Row Tx %d 失败 - %s\n', tx, ME.message);
        end
        
        if mod(tx, 4) == 0
            fprintf('    Row Tx: %d/%d\n', tx, N);
        end
    end
    
    %% Column发射 -> Row接收
    rf_CR = cell(N, 1);
    tstart_CR = zeros(N, 1);
    
    for tx = 1:N
        apo_tx = zeros(1, N);
        apo_tx(tx) = 1;
        xdc_apodization(Tx_col, 0, apo_tx);
        xdc_focus(Tx_col, 0, [0 0 500e-3]);
        
        xdc_apodization(Rx_row, 0, ones(1, N));
        
        try
            [rf, ts] = calc_scat_batched(Tx_col, Rx_row, positions, amplitudes, params.fs);
            
            if ~isempty(rf) && numel(rf) > 50
                rf_CR{tx} = single(rf);
                tstart_CR(tx) = ts;
            end
        catch ME
            fprintf('    警告: Col Tx %d 失败 - %s\n', tx, ME.message);
        end
        
        if mod(tx, 4) == 0
            fprintf('    Col Tx: %d/%d\n', tx, N);
        end
    end
    
    fprintf('    有效数据: RC=%d, CR=%d\n', sum(~cellfun(@isempty, rf_RC)), ...
            sum(~cellfun(@isempty, rf_CR)));
end

function [rf_total, tstart_total] = calc_scat_batched(Tx, Rx, pos, amp, fs)
%CALC_SCAT_BATCHED 分批计算散射以避免Field II内存溢出
    
    BATCH_SIZE = 1000; % 关键参数：每批次处理的点数，显存不足可继续调小
    n_points = size(pos, 1);
    
    % 如果点数很少，直接计算，无需分批
    if n_points <= BATCH_SIZE
        [rf_total, tstart_total] = calc_scat(Tx, Rx, pos, amp);
        return;
    end
    
    % 初始化分批变量
    n_batches = ceil(n_points / BATCH_SIZE);
    rf_parts = cell(n_batches, 1);
    tstart_parts = zeros(n_batches, 1);
    valid_mask = false(n_batches, 1);
    
    % 循环计算每一批
    for b = 1:n_batches
        idx_start = (b-1)*BATCH_SIZE + 1;
        idx_end = min(b*BATCH_SIZE, n_points);
        
        pos_batch = pos(idx_start:idx_end, :);
        amp_batch = amp(idx_start:idx_end);
        
        try
            [rf, ts] = calc_scat(Tx, Rx, pos_batch, amp_batch);
            if ~isempty(rf)
                rf_parts{b} = rf;
                tstart_parts(b) = ts;
                valid_mask(b) = true;
            end
        catch
            % 忽略空批次或失败批次
        end
    end
    
    if ~any(valid_mask)
        rf_total = [];
        tstart_total = 0;
        return;
    end
    
    % --- 合并结果（核心逻辑：按时间对齐叠加）---
    valid_indices = find(valid_mask);
    
    % 1. 找到最早的起始时间
    min_tstart = min(tstart_parts(valid_indices));
    tstart_total = min_tstart;
    
    % 2. 计算合成信号所需的总长度
    max_sample_index = 0;
    for i = 1:length(valid_indices)
        idx = valid_indices(i);
        % 计算当前批次相对于最早时间的样本偏移量
        start_sample_offset = round((tstart_parts(idx) - min_tstart) * fs);
        end_sample_index = start_sample_offset + size(rf_parts{idx}, 1);
        if end_sample_index > max_sample_index
            max_sample_index = end_sample_index;
        end
    end
    
    % 3. 初始化总RF数据并叠加
    n_channels = size(rf_parts{valid_indices(1)}, 2);
    rf_total = zeros(max_sample_index, n_channels);
    
    for i = 1:length(valid_indices)
        idx = valid_indices(i);
        rf = rf_parts{idx};
        start_sample_offset = round((tstart_parts(idx) - min_tstart) * fs);
        
        s_idx = start_sample_offset + 1;
        e_idx = start_sample_offset + size(rf, 1);
        
        rf_total(s_idx:e_idx, :) = rf_total(s_idx:e_idx, :) + rf;
    end
end
