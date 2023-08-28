function y_p = Phi_p_times_x_p(Phi_p, N, x_p)
% Input: periodic part Phi_p, row dimension of Phi N, periodic part x_p
% The function compute Phi*x using period and Toeplitz properties
% Output: periodic part of y: y_p = y(1:p)
    
    p = size(Phi_p, 1);
    k = floor(N / p);
    r = mod(N, p);
    y_p = k * FFT_Toep_matrix_multiplication(Phi_p, x_p);
    y_p = y_p + FFT_Toep_matrix_multiplication(Phi_p(:,1:r), x_p(1:r));

end

