function [all_positions, all_amplitudes] = update_blood_positions(phantom, params, t_current)
%UPDATE_BLOOD_POSITIONS 更新血液散射点位置
%   实现Poiseuille抛物线流速分布

    % 归一化半径
    normalized_r = phantom.r_blood / params.vessel_radius;
    
    % Poiseuille流速分布: v(r) = v_max * (1 - (r/R)^2)
    v_local = params.v_max * (1 - normalized_r.^2);
    
    % Y方向位移
    y_blood = phantom.y_blood_init + v_local * t_current;
    
    % 周期性边界条件
    y_blood = mod(y_blood + params.vessel_length/2, params.vessel_length) - params.vessel_length/2;
    
    % 计算血液3D位置
    blood_positions = zeros(phantom.n_blood, 3);
    blood_positions(:,1) = phantom.r_blood .* cos(phantom.theta_blood);
    blood_positions(:,2) = y_blood;
    blood_positions(:,3) = phantom.r_blood .* sin(phantom.theta_blood) + params.vessel_depth;
    
    % 合并静态和动态散射点
    all_positions = [phantom.static_positions; blood_positions];
    all_amplitudes = [phantom.static_amplitudes; phantom.blood_amplitudes];
end