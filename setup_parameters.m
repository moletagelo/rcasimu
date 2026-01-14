function params = setup_parameters()
%SETUP_PARAMETERS 设置所有仿真参数
%   减少参数值以降低内存使用

    %% 物理参数
    params.c = 1540;                    % 声速 (m/s)
    params.f0 = 10e6;                    % 中心频率 (Hz)
    params.fs = 100e6;                  % 采样频率 (Hz)
    params.lambda = params.c / params.f0;
    
    %% 阵列参数 (减少阵元数以节省内存)
    params.N_elements = 32;             % 每个方向阵元数 
    params.pitch = params.lambda;       % 阵元间距
    params.kerf = 0.01e-3;             % 阵元间隙
    params.element_size = params.pitch - params.kerf;
    params.array_size = params.N_elements * params.pitch;
    
    % 阵元位置
    params.row_pos_x = ((1:params.N_elements) - (params.N_elements+1)/2) * params.pitch;
    params.col_pos_y = ((1:params.N_elements) - (params.N_elements+1)/2) * params.pitch;
    
    %% 血管参数
    params.vessel_radius = 2.5e-3;      % 血管半径 2.5mm
    params.vessel_length = 15e-3;       % 血管长度 15mm
    params.vessel_depth = 30e-3;        % 血管深度 30mm
    
    %% 血流参数
    params.v_max = 0.4;                 % 最大流速 m/s
    params.PRF = 2000;                  % 脉冲重复频率
    params.dt_frame = 1/params.PRF;
    
    %% 散射点数量 (大幅减少以节省内存)
    params.n_blood = 500;               % 血液散射点 
    params.n_wall = 300;                % 血管壁散射点 
    params.n_tissue = 800;              % 组织散射点 
    
    %% 仿真帧数
    params.n_frames = 10;                % 帧数 
    
    %% 成像网格 (减少分辨率以节省内存)
    params.n_x = 32;                    % X方向像素 (原48改为32)
    params.n_y = 32;                    % Y方向像素
    params.n_z = 40;                    % Z方向像素 (原64改为40)
    
    params.x_range = linspace(-10e-3, 10e-3, params.n_x);
    params.y_range = linspace(-10e-3, 10e-3, params.n_y);
    params.z_range = linspace(22e-3, 40e-3, params.n_z);
    
    %% 显示参数
    fprintf('  阵列: %dx%d 阵元\n', params.N_elements, params.N_elements);
    fprintf('  血管: 半径%.1fmm, 深度%.0fmm\n', params.vessel_radius*1e3, params.vessel_depth*1e3);
    fprintf('  散射点: 壁%d + 组织%d + 血液%d\n', params.n_wall, params.n_tissue, params.n_blood);
    fprintf('  成像网格: %dx%dx%d\n', params.n_x, params.n_y, params.n_z);
end