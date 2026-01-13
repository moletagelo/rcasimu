function [rf_RC, tstart_RC, rf_CR, tstart_CR] = acquire_rf_data(...
    Tx_row, Rx_col, Tx_col, Rx_row, positions, amplitudes, params)
%ACQUIRE_RF_DATA 采集Row-Column RF数据
%   分批次采集以节省内存

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
        
        % 计算散射响应
        try
            [rf, ts] = calc_scat(Tx_row, Rx_col, positions, amplitudes);
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
            [rf, ts] = calc_scat(Tx_col, Rx_row, positions, amplitudes);
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