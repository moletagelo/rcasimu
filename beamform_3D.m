function volume = beamform_3D(rf_RC, tstart_RC, rf_CR, tstart_CR, params)
%BEAMFORM_3D 三维延迟叠加波束形成
%   使用向量化操作加速

    N = params.N_elements;
    c = params.c;
    fs = params.fs;
    
    n_x = params.n_x;
    n_y = params.n_y;
    n_z = params.n_z;
    
    % 输出体积 (使用单精度)
    volume = zeros(n_y, n_x, n_z, 'single');
    
    % 创建3D网格
    [X, Y, Z] = ndgrid(single(params.x_range), single(params.y_range), single(params.z_range));
    
    %% Row->Column 波束形成
    for tx = 1:N
        if isempty(rf_RC{tx})
            continue;
        end
        
        rf_raw = rf_RC{tx};
        t0 = tstart_RC(tx);
        
        % 确保2D
        if isvector(rf_raw)
            rf_raw = rf_raw(:);
        end
        
        [n_samp, n_ch] = size(rf_raw);
        if n_samp < 50
            continue;
        end
        
        % 包络检测
        rf_env = single(abs(hilbert(double(rf_raw))));
        
        % Row发射阵元X位置
        tx_x = params.row_pos_x(tx);
        
        % 发射距离 (Row阵元沿Y延伸，无Y聚焦)
        tx_dist = sqrt((X - tx_x).^2 + Z.^2);
        
        for rx = 1:min(n_ch, N)
            % Column接收阵元Y位置
            rx_y = params.col_pos_y(rx);
            
            % 接收距离 (Column阵元沿X延伸，无X聚焦)
            rx_dist = sqrt((Y - rx_y).^2 + Z.^2);
            
            % 总飞行时间
            tof = (tx_dist + rx_dist) / c;
            sample_idx = round((tof - t0) * fs) + 1;
            
            % 有效范围
            valid = (sample_idx >= 1) & (sample_idx <= n_samp);
            sample_idx(~valid) = 1;
            
            % 提取通道数据
            rf_ch = rf_env(:, rx);
            
            % 插值取值并累加
            values = rf_ch(sample_idx);
            values(~valid) = 0;
            
            volume = volume + permute(values, [2 1 3]);
        end
    end
    
    %% Column->Row 波束形成
    for tx = 1:N
        if isempty(rf_CR{tx})
            continue;
        end
        
        rf_raw = rf_CR{tx};
        t0 = tstart_CR(tx);
        
        if isvector(rf_raw)
            rf_raw = rf_raw(:);
        end
        
        [n_samp, n_ch] = size(rf_raw);
        if n_samp < 50
            continue;
        end
        
        rf_env = single(abs(hilbert(double(rf_raw))));
        
        % Column发射阵元Y位置
        tx_y = params.col_pos_y(tx);
        
        % 发射距离
        tx_dist = sqrt((Y - tx_y).^2 + Z.^2);
        
        for rx = 1:min(n_ch, N)
            % Row接收阵元X位置
            rx_x = params.row_pos_x(rx);
            
            % 接收距离
            rx_dist = sqrt((X - rx_x).^2 + Z.^2);
            
            tof = (tx_dist + rx_dist) / c;
            sample_idx = round((tof - t0) * fs) + 1;
            
            valid = (sample_idx >= 1) & (sample_idx <= n_samp);
            sample_idx(~valid) = 1;
            
            rf_ch = rf_env(:, rx);
            values = rf_ch(sample_idx);
            values(~valid) = 0;
            
            volume = volume + permute(values, [2 1 3]);
        end
    end
    
    volume = double(volume);
end