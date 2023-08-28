function y_p = Phi_p_times_x(Phi_p, N, x)
% Input: periodic part Phi_p, row dimension of Phi N, N-dim vector x
% The function compute Phi'*x using period and Toeplitz properties
% Output: periodic part of y: y_p = y(1:p)
    
    p = size(Phi_p, 1);
    k = floor(N / p);
    r = mod(N, p);
    y_p = zeros(p, 1);
    
    for i = 1:k
        y_p = y_p + FFT_Toep_matrix_multiplication(Phi_p, x((i-1)*p+1:i*p));
    end
    if r ~= 0
    y_p = y_p + FFT_Toep_matrix_multiplication(Phi_p(:,1:r), x(N-r+1:N));
    end
end

