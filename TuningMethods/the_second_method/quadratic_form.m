function res = quadratic_form(Phi_p, N, n, x, Ut, Vt, alpha)
% Input: period part Phi_p, dimension of Phi: N*n, period p,
%        N-dim x, semi-separable representation U & V, and alpha
% H = Phi * K * Phi' + alpha * I
% Phi is periodic Toeplitz matrix with period p. K = S(U,V) with rank k.
% Output: x'* H * x

    % y = Phi' * x
    y_p = Phi_p_times_x(Phi_p', N, x);
    
    % y'*K*y
    y = period2whole_vector(y_p, n);
    res1 = y' * semi_vector_multiplication(Ut, Vt, y);
    
    % alpha*x'*x
    res2 = alpha * (x'*x);
    
    % Output result
    res = res1 + res2;
end

