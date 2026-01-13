function visualize_results(volume_data, flow_volume, params)
%VISUALIZE_RESULTS 可视化三维成像结果

    % 选择帧
    frame = params.n_frames;
    vol = volume_data(:,:,:,frame);
    flow = flow_volume(:,:,:,frame);
    
    % 归一化
    vol = vol / (max(vol(:)) + eps);
    flow = flow / (max(flow(:)) + eps);
    
    % dB转换
    vol_dB = 20*log10(vol + 1e-6);
    vol_dB = max(vol_dB, -50);
    
    % 切片索引
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
    colormap(gca, gray(256)); caxis([-50 0]);
    colorbar; xlabel('X (mm)'); ylabel('Z (mm)');
    title('X-Z切片 (血管横截面)');
    axis image; set(gca, 'YDir', 'reverse');
    
    % Y-Z切片
    subplot(1,3,2);
    img = squeeze(vol_dB(:,ix_mid,:))';
    imagesc(y_mm, z_mm, img);
    hold on;
    plot([-1 1]*params.vessel_length/2*1e3, [1 1]*(params.vessel_depth - params.vessel_radius)*1e3, 'r--', 'LineWidth', 1.5);
    plot([-1 1]*params.vessel_length/2*1e3, [1 1]*(params.vessel_depth + params.vessel_radius)*1e3, 'r--', 'LineWidth', 1.5);
    hold off;
    colormap(gca, gray(256)); caxis([-50 0]);
    colorbar; xlabel('Y (mm)'); ylabel('Z (mm)');
    title('Y-Z切片 (血管纵向)');
    axis image; set(gca, 'YDir', 'reverse');
    
    % X-Y切片
    subplot(1,3,3);
    img = squeeze(vol_dB(:,:,iz_vessel));
    imagesc(x_mm, y_mm, img);
    hold on;
    viscircles([0 0], params.vessel_radius*1e3, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    hold off;
    colormap(gca, gray(256)); caxis([-50 0]);
    colorbar; xlabel('X (mm)'); ylabel('Y (mm)');
    title(sprintf('X-Y切片 (Z=%.0fmm)', params.vessel_depth*1e3));
    axis image;
    
    sgtitle('Row-Column三维血管成像结果', 'FontSize', 14, 'FontWeight', 'bold');
    
    %% 图2: MIP投影
    figure('Name', 'MIP投影', 'Position', [100 80 1200 400], 'Color', 'w');
    
    subplot(1,3,1);
    mip = squeeze(max(vol_dB, [], 1));
    imagesc(y_mm, z_mm, mip');
    colormap(gca, gray(256)); caxis([-50 0]); colorbar;
    xlabel('Y (mm)'); ylabel('Z (mm)');
    title('MIP (沿X)'); axis image;
    set(gca, 'YDir', 'reverse');
    
    subplot(1,3,2);
    mip = squeeze(max(vol_dB, [], 2));
    imagesc(x_mm, z_mm, mip');
    colormap(gca, gray(256)); caxis([-50 0]); colorbar;
    xlabel('X (mm)'); ylabel('Z (mm)');
    title('MIP (沿Y)'); axis image;
    set(gca, 'YDir', 'reverse');
    
    subplot(1,3,3);
    mip = squeeze(max(vol_dB, [], 3));
    imagesc(x_mm, y_mm, mip);
    colormap(gca, gray(256)); caxis([-50 0]); colorbar;
    xlabel('X (mm)'); ylabel('Y (mm)');
    title('MIP (沿Z, C-scan)'); axis image;
    
    sgtitle('最大强度投影', 'FontSize', 14);
    
    %% 图3: B模式与血流
    figure('Name', 'B模式与血流', 'Position', [150 60 1200 500], 'Color', 'w');
    
    % B模式
    subplot(2,3,1);
    imagesc(x_mm, z_mm, squeeze(vol_dB(iy_mid,:,:))');
    colormap(gca, gray); caxis([-50 0]);
    xlabel('X'); ylabel('Z'); title('B模式 X-Z');
    axis image; set(gca, 'YDir', 'reverse');
    
    subplot(2,3,2);
    imagesc(y_mm, z_mm, squeeze(vol_dB(:,ix_mid,:))');
    colormap(gca, gray); caxis([-50 0]);
    xlabel('Y'); ylabel('Z'); title('B模式 Y-Z');
    axis image; set(gca, 'YDir', 'reverse');
    
    subplot(2,3,3);
    imagesc(x_mm, y_mm, squeeze(vol_dB(:,:,iz_vessel)));
    colormap(gca, gray); caxis([-50 0]);
    xlabel('X'); ylabel('Y'); title('B模式 X-Y');
    axis image;
    
    % 血流
    subplot(2,3,4);
    imagesc(x_mm, z_mm, squeeze(flow(iy_mid,:,:))');
    colormap(gca, hot);
    xlabel('X'); ylabel('Z'); title('血流 X-Z');
    axis image; set(gca, 'YDir', 'reverse');
    
    subplot(2,3,5);
    imagesc(y_mm, z_mm, squeeze(flow(:,ix_mid,:))');
    colormap(gca, hot);
    xlabel('Y'); ylabel('Z'); title('血流 Y-Z');
    axis image; set(gca, 'YDir', 'reverse');
    
    subplot(2,3,6);
    imagesc(x_mm, y_mm, squeeze(flow(:,:,iz_vessel)));
    colormap(gca, hot);
    xlabel('X'); ylabel('Y'); title('血流 X-Y');
    axis image;
    
    sgtitle('B模式与血流对比', 'FontSize', 14);
    
    %% 图4: 3D可视化
    figure('Name', '三维可视化', 'Position', [200 50 900 700], 'Color', 'w');
    
    [Xm, Ym, Zm] = meshgrid(x_mm, y_mm, z_mm);
    
    subplot(2,2,1);
    slice(Xm, Ym, Zm, vol, 0, 0, params.vessel_depth*1e3);
    shading interp; colormap(gca, gray);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D切片'); axis equal tight;
    set(gca, 'ZDir', 'reverse'); view(30, 25);
    
    subplot(2,2,2);
    thresh = 0.2;
    try
        fv = isosurface(Xm, Ym, Zm, vol, thresh);
        if ~isempty(fv.vertices)
            p = patch(fv);
            set(p, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            lighting gouraud; camlight;
        end
    catch
    end
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('等值面'); axis equal tight;
    set(gca, 'ZDir', 'reverse'); view(30, 25); grid on;
    
    subplot(2,2,3);
    slice(Xm, Ym, Zm, flow, 0, 0, params.vessel_depth*1e3);
    shading interp; colormap(gca, hot);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('血流3D'); axis equal tight;
    set(gca, 'ZDir', 'reverse'); view(30, 25);
    
    subplot(2,2,4);
    % 剖面分析
    profile_x = squeeze(vol(iy_mid, :, iz_vessel));
    profile_x = profile_x / max(profile_x);
    plot(x_mm, profile_x, 'b-', 'LineWidth', 2);
    hold on;
    xline(-params.vessel_radius*1e3, 'r--');
    xline(params.vessel_radius*1e3, 'r--');
    hold off;
    xlabel('X (mm)'); ylabel('强度');
    title('X方向剖面'); grid on;
    legend('成像', '血管边界');
    
    sgtitle('三维可视化', 'FontSize', 14);
    
    fprintf('可视化完成\n');
end