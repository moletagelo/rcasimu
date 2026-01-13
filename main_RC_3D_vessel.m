%% =========================================================================
%  Row-Column阵列三维动态血管成像 - 主程序
%  优化版：减少内存使用，拆分模块
%  =========================================================================
clear all;
close all;
clc;

fprintf('==============================================\n');
fprintf('  Row-Column 三维动态血管成像仿真\n');
fprintf('==============================================\n\n');

%% ========== 步骤1: 设置参数 ==========
fprintf('[1/6] 设置参数...\n');
params = setup_parameters();

%% ========== 步骤2: 初始化Field II ==========
fprintf('[2/6] 初始化Field II...\n');
field_init(-1);
set_field('fs', params.fs);
set_field('c', params.c);
fprintf('  Field II初始化完成\n');

%% ========== 步骤3: 创建换能器 ==========
fprintf('[3/6] 创建换能器...\n');
[Tx_row, Rx_col, Tx_col, Rx_row] = create_transducers(params);
fprintf('  换能器创建完成\n');

%% ========== 步骤4: 创建散射点模型 ==========
fprintf('[4/6] 创建血管散射点模型...\n');
phantom = create_phantom(params);
fprintf('  散射点: %d个\n', phantom.n_total);

%% ========== 步骤5: 多帧采集与重建 ==========
fprintf('[5/6] 开始多帧采集与重建...\n\n');

% 初始化存储
volume_data = zeros(params.n_y, params.n_x, params.n_z, params.n_frames);
flow_volume = zeros(params.n_y, params.n_x, params.n_z, params.n_frames);

for frame = 1:params.n_frames
    fprintf('===== 帧 %d/%d =====\n', frame, params.n_frames);
    frame_tic = tic;
    
    % 计算当前时刻血液位置
    t_current = (frame - 1) * params.dt_frame;
    [all_positions, all_amplitudes] = update_blood_positions(phantom, params, t_current);
    
    % 采集RF数据
    fprintf('  采集RF数据...\n');
    [rf_RC, tstart_RC, rf_CR, tstart_CR] = acquire_rf_data(...
        Tx_row, Rx_col, Tx_col, Rx_row, all_positions, all_amplitudes, params);
    
    % 三维波束形成
    fprintf('  波束形成...\n');
    vol_frame = beamform_3D(rf_RC, tstart_RC, rf_CR, tstart_CR, params);
    
    volume_data(:,:,:,frame) = vol_frame;
    
    % 计算血流(帧间差分)
    if frame > 1
        flow_volume(:,:,:,frame) = abs(vol_frame - volume_data(:,:,:,frame-1));
    end
    
    fprintf('  帧 %d 完成, 耗时 %.1f 秒\n\n', frame, toc(frame_tic));
    
    % 清理临时变量释放内存
    clear rf_RC rf_CR all_positions all_amplitudes;
end

% 第一帧血流
flow_volume(:,:,:,1) = flow_volume(:,:,:,min(2, params.n_frames));

%% ========== 步骤6: 可视化 ==========
fprintf('[6/6] 生成可视化...\n');
visualize_results(volume_data, flow_volume, params);

%% ========== 清理 ==========
fprintf('\n清理资源...\n');
xdc_free(Tx_row);
xdc_free(Rx_col);
xdc_free(Tx_col);
xdc_free(Rx_row);
field_end;

%% ========== 保存结果 ==========
fprintf('保存结果...\n');
save('RC_3D_results.mat', 'volume_data', 'flow_volume', 'params', '-v7.3');

fprintf('\n========== 仿真完成 ==========\n');