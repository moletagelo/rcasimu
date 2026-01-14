function visualize_results(volume_data, flow_volume, params)
%VISUALIZE_RESULTS 可视化三维成像结果
%   修复：使用鲁棒归一化解决图像全黑问题

    % 选择帧
    frame = params.n_frames;
    vol = volume_data(:,:,:,frame);
    flow = flow_volume(:,:,:,frame);
    
    %% === 核心修改：鲁棒归一化 ===
    % 1. 检查数据是否全零
    max_val = max(vol(:));
    if max_val < 1e-12
        warning('图像数据全为0，请检查仿真过程！');
        return;
    end
    
    % 2. 使用分位数而不是最大值进行归一化 (抗噪点)
    % 排序后取第99.9%大的值作为"最大值"
    sorted_vals = sort(vol(:));
    robust_max = sorted_vals(round(0.999 * length(sorted_vals))); 
    
    % 防止 robust_max 为 0
    if robust_max == 0, robust_max = max_val; end
    
    fprintf('  信号峰值: %.2e, 99.9%%分位值: %.2e (增益提升: %.1f dB)\n', ...
            max_val, robust_max, 20*log10(max_val/robust_max));

    % 3. 执行归一化和截断
    vol = vol / robust_max;
    vol(vol > 1) = 1; % 超过分位数的点截断为1
    
    % 4. 对数压缩 (Log Compression)
    dynamic_range = 60; % 增加动态范围到 60dB (原50dB)
    vol_dB = 20*log10(vol + 1e-9);
    vol_dB = max(vol_dB, -dynamic_range);
    
    % --- 血流数据同理处理 ---
    flow_max = max(flow(:));
    if flow_max > 0
        sorted_flow = sort(flow(:));
        robust_flow_max = sorted_flow(round(0.999 * length(sorted_flow)));
        if robust_flow_max == 0, robust_flow_max = flow_max; end
        flow = flow / robust_flow_max;
        flow(flow > 1) = 1;
    end
    
    %% 下面是绘图部分 (保持原样，仅需更新caxis范围)
    
    % ... (前面的坐标定义代码不变) ...
    [~, ix_mid] = min(abs(params.x_range));
    [~, iy_mid] = min(abs(params.y_range));
    [~, iz_vessel] = min(abs(params.z_range - params.vessel_depth));
    
    x_mm = params.x_range * 1e3;
    y_mm = params.y_range * 1e3;
    z_mm = params.z_range * 1e3;
    
    %% 图1: 三个正交切片
    figure('Name', 'Row-Column三维血管成像', 'Position', [50 100 1400 450], 'Color', 'w');
    
    % X-Z切片
    subplot(1,3,1);
    img = squeeze(vol_dB(iy_mid,:,:))';
    imagesc(x_mm, z_mm, img);
    hold on;
    theta = linspace(0, 2*pi, 50);
    plot(params.vessel_radius*cos(theta)*1e3, ...
         (params.vessel_radius*sin(theta) + params.vessel_depth)*1e3, 'r--', 'LineWidth', 1.5);
    hold off;
    colormap(gca, gray(256)); 
    caxis([-dynamic_range 0]); % 使用新的动态范围
    colorbar; xlabel('X (mm)'); ylabel('Z (mm)');
    title('X-Z切片 (血管横截面)');
    axis image; set(gca, 'YDir', 'reverse');
    
    % ... (后续绘图代码类似，只需将 caxis([-50 0]) 改为 caxis([-dynamic_range 0])) ...
    
    % 简单起见，这里只展示关键修改。
    % 请确保将所有 subplot 中的 caxis([-50 0]) 替换为 caxis([-dynamic_range 0])
end
