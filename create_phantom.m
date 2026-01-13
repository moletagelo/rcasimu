function phantom = create_phantom(params)
%CREATE_PHANTOM 创建血管散射点模型
%   返回结构体包含静态散射点和血液初始位置

    rng(123);  % 固定随机种子
    
    %% 血管壁散射点 (圆柱面)
    n_wall = params.n_wall;
    theta_wall = 2*pi*rand(n_wall, 1);
    y_wall = params.vessel_length * (rand(n_wall, 1) - 0.5);
    
    wall_positions = zeros(n_wall, 3);
    wall_positions(:,1) = params.vessel_radius * cos(theta_wall);
    wall_positions(:,2) = y_wall;
    wall_positions(:,3) = params.vessel_radius * sin(theta_wall) + params.vessel_depth;
    wall_amplitudes = 5 * ones(n_wall, 1);
    
    %% 周围组织散射点
    n_tissue = params.n_tissue;
    tissue_range_xy = 12e-3;
    tissue_z_min = 22e-3;
    tissue_z_max = 40e-3;
    
    tissue_x = tissue_range_xy * (rand(n_tissue, 1) - 0.5);
    tissue_y = tissue_range_xy * (rand(n_tissue, 1) - 0.5);
    tissue_z = (tissue_z_max - tissue_z_min) * rand(n_tissue, 1) + tissue_z_min;
    
    % 移除血管内部点
    dist_from_axis = sqrt(tissue_x.^2 + (tissue_z - params.vessel_depth).^2);
    in_vessel_y = abs(tissue_y) < params.vessel_length/2;
    inside = (dist_from_axis < params.vessel_radius * 1.3) & in_vessel_y;
    
    tissue_x(inside) = [];
    tissue_y(inside) = [];
    tissue_z(inside) = [];
    
    tissue_positions = [tissue_x, tissue_y, tissue_z];
    tissue_amplitudes = 0.5 * ones(size(tissue_positions, 1), 1);
    
    %% 合并静态散射点
    phantom.static_positions = [wall_positions; tissue_positions];
    phantom.static_amplitudes = [wall_amplitudes; tissue_amplitudes];
    
    %% 血液初始化
    n_blood = params.n_blood;
    phantom.theta_blood = 2*pi*rand(n_blood, 1);
    phantom.r_blood = params.vessel_radius * sqrt(rand(n_blood, 1)) * 0.88;
    phantom.y_blood_init = params.vessel_length * (rand(n_blood, 1) - 0.5);
    phantom.blood_amplitudes = 0.15 * ones(n_blood, 1);
    
    phantom.n_static = size(phantom.static_positions, 1);
    phantom.n_blood = n_blood;
    phantom.n_total = phantom.n_static + n_blood;
    
    fprintf('  壁散射点: %d\n', n_wall);
    fprintf('  组织散射点: %d\n', size(tissue_positions, 1));
    fprintf('  血液散射点: %d\n', n_blood);
end

function [all_positions, all_amplitudes] = update_blood_positions(phantom, params, t_current)
%UPDATE_BLOOD_POSITIONS 更新血液散射点位置
    
    % Poiseuille流速分布
    normalized_r = phantom.r_blood / params.vessel_radius;
    v_local = params.v_max * (1 - normalized_r.^2);
    
    % Y方向位移
    y_blood = phantom.y_blood_init + v_local * t_current;
    y_blood = mod(y_blood + params.vessel_length/2, params.vessel_length) - params.vessel_length/2;
    
    % 血液位置
    blood_positions = zeros(phantom.n_blood, 3);
    blood_positions(:,1) = phantom.r_blood .* cos(phantom.theta_blood);
    blood_positions(:,2) = y_blood;
    blood_positions(:,3) = phantom.r_blood .* sin(phantom.theta_blood) + params.vessel_depth;
    
    % 合并
    all_positions = [phantom.static_positions; blood_positions];
    all_amplitudes = [phantom.static_amplitudes; phantom.blood_amplitudes];
end