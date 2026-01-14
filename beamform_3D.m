function volume = beamform_3D(rf_RC, tstart_RC, rf_CR, tstart_CR, params)
%BEAMFORM_3D_GPU 使用GPU加速的三维波束形成
%   利用 MATLAB gpuArray 加速大规模矩阵运算

    % 检查GPU是否可用
    try
        g = gpuDevice;
        use_gpu = true;
        fprintf('  [GPU加速] 使用设备: %s\n', g.Name);
    catch
        use_gpu = false;
        fprintf('  [CPU模式] 未检测到GPU，回退到CPU\n');
    end

    N = params.N_elements;
    c = params.c;
    fs = params.fs;
    
    % 将坐标网格移入 GPU
    x_range = single(params.x_range);
    y_range = single(params.y_range);
    z_range = single(params.z_range);
    
    if use_gpu
        x_range = gpuArray(x_range);
        y_range = gpuArray(y_range);
        z_range = gpuArray(z_range);
        volume = gpuArray.zeros(params.n_y, params.n_x, params.n_z, 'single');
    else
        volume = zeros(params.n_y, params.n_x, params.n_z, 'single');
    end
    
    [X, Y, Z] = ndgrid(x_range, y_range, z_range);
    
    %% Row -> Column (发射X位置固定，接收Y位置固定)
    row_pos_x = params.row_pos_x; 
    col_pos_y = params.col_pos_y;
    
    if use_gpu
        row_pos_x = gpuArray(single(row_pos_x));
        col_pos_y = gpuArray(single(col_pos_y));
    end

    fprintf('  正在进行 Row->Col GPU波束形成...\n');
    
    for tx = 1:N
        if isempty(rf_RC{tx}), continue; end
        
        % 数据预处理 (CPU端)
        rf_raw = rf_RC{tx};
        t0 = tstart_RC(tx);
        
        % 包络检测建议在CPU做，或者移到GPU
        rf_env = single(abs(hilbert(double(rf_raw))));
        
        if use_gpu
            rf_env = gpuArray(rf_env);
        end
        
        [n_samp, n_ch] = size(rf_env);
        
        % --- GPU 矢量化计算 ---
        % 发射距离场 (仅与X,Z有关)
        tx_x = row_pos_x(tx);
        Tx_Dist = sqrt((X - tx_x).^2 + Z.^2);
        
        % 并行处理所有接收通道 (显存允许的情况下)
        % 为了防止显存爆炸，我们还是逐通道循环，但在GPU上计算是大矩阵操作
        
        for rx = 1:min(n_ch, N)
            rx_y = col_pos_y(rx);
            
            % 接收距离场 (仅与Y,Z有关)
            Rx_Dist = sqrt((Y - rx_y).^2 + Z.^2);
            
            % 延迟计算
            TOF = (Tx_Dist + Rx_Dist) / c;
            Sample_Idx = round((TOF - t0) * fs) + 1;
            
            % 边界检查
            Mask = (Sample_Idx >= 1) & (Sample_Idx <= n_samp);
            Sample_Idx(~Mask) = 1;
            
            % 索引取值
            RF_Data = rf_env(:, rx);
            Val = RF_Data(Sample_Idx);
            Val(~Mask) = 0;
            
            % 累加
            volume = volume + Val;
        end
    end
    
    %% Column -> Row
    fprintf('  正在进行 Col->Row GPU波束形成...\n');
    
    for tx = 1:N
        if isempty(rf_CR{tx}), continue; end
        
        rf_raw = rf_CR{tx};
        t0 = tstart_CR(tx);
        rf_env = single(abs(hilbert(double(rf_raw))));
        
        if use_gpu, rf_env = gpuArray(rf_env); end
        [n_samp, n_ch] = size(rf_env);
        
        tx_y = col_pos_y(tx);
        Tx_Dist = sqrt((Y - tx_y).^2 + Z.^2);
        
        for rx = 1:min(n_ch, N)
            rx_x = row_pos_x(rx);
            Rx_Dist = sqrt((X - rx_x).^2 + Z.^2);
            
            TOF = (Tx_Dist + Rx_Dist) / c;
            Sample_Idx = round((TOF - t0) * fs) + 1;
            
            Mask = (Sample_Idx >= 1) & (Sample_Idx <= n_samp);
            Sample_Idx(~Mask) = 1;
            
            RF_Data = rf_env(:, rx);
            Val = RF_Data(Sample_Idx);
            Val(~Mask) = 0;
            
            volume = volume + Val;
        end
    end
    
    % 从GPU取回数据
    if use_gpu
        volume = gather(volume);
    end
    
    volume = double(volume);
end
