function [Tx_row, Rx_col, Tx_col, Rx_row] = create_transducers(params)
%CREATE_TRANSDUCERS 创建Row-Column换能器阵列
%   Tx_row: Row发射阵列
%   Rx_col: Column接收阵列
%   Tx_col: Column发射阵列
%   Rx_row: Row接收阵列

    N = params.N_elements;
    
    %% 脉冲响应和激励
    t_ir = (-1/params.f0 : 1/params.fs : 1/params.f0);
    impulse_response = sin(2*pi*params.f0*t_ir) .* hanning(length(t_ir))';
    
    t_ex = (0 : 1/params.fs : 2/params.f0);
    excitation = sin(2*pi*params.f0*t_ex);
    
    %% Row阵列 (沿X方向排列，每个阵元沿Y延伸)
    % 使用xdc_linear_array
    Tx_row = xdc_linear_array(N, params.element_size, params.array_size, ...
                               params.kerf, 1, 4, [0 0 50e-3]);
    xdc_impulse(Tx_row, impulse_response);
    xdc_excitation(Tx_row, excitation);
    
    Rx_row = xdc_linear_array(N, params.element_size, params.array_size, ...
                               params.kerf, 1, 4, [0 0 50e-3]);
    xdc_impulse(Rx_row, impulse_response);
    
    %% Column阵列 (沿Y方向排列，每个阵元沿X延伸)
    % 手动创建矩形阵元
    rect_col = zeros(N, 19);
    for i = 1:N
        rect_col(i, 1) = i;                             % 阵元编号
        rect_col(i, 2) = params.array_size;             % X尺寸(宽)
        rect_col(i, 3) = params.element_size;           % Y尺寸(窄)
        rect_col(i, 4:6) = [4 1 1];                     % 子阵元划分
        rect_col(i, 7) = 0;                             % X中心
        rect_col(i, 8) = params.col_pos_y(i);           % Y中心
        rect_col(i, 9) = 0;                             % Z位置
        rect_col(i, 11:19) = [0 0 1 1 0 0 0 0 0];
    end
    
    Rx_col = xdc_rectangles(rect_col, [0 0 50e-3], [0 0 1e9]);
    xdc_impulse(Rx_col, impulse_response);
    
    Tx_col = xdc_rectangles(rect_col, [0 0 50e-3], [0 0 1e9]);
    xdc_impulse(Tx_col, impulse_response);
    xdc_excitation(Tx_col, excitation);
    
    fprintf('  Row阵列: %d阵元, Column阵列: %d阵元\n', N, N);
end