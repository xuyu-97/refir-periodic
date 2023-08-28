function res = method1_quadratic_form(Phi_p, N, n, x, Ut, Vt, alpha)
%  Input: period part Phi_p, dimension of Phi: N*n, period p,
%        N-dim x, semi-separable representation U & V, and alpha
% H = Phi * K * Phi' + alpha * I
% Phi is periodic Toeplitz matrix with period p. K = S(U,V) with rank k.
% Output: x'* H * x
% Method: Phi = QR (full), x'*H*x = x'*Q*(R*K*R'+alpha*I)*Q'*x
    
    p = size(Phi_p, 1);
    [Qp, Rp] = HouseholderQR_periodic(Phi_p, N);
    
    % Q'*x = Pp*...*P1*x
    for k = 1:p
        b_k = Qp(p+2,k);
        v_k = Qp(1:(p+1),k);
        v_k = [v_k(1); period2whole_vector(v_k(2:(p+1)), N-k)];
        x(k:N) = x(k:N) - b_k*v_k*(v_k'*x(k:N));
    end
    
    % y = R'*x, y is periodic
    y = Rp'*x(1:p);
    y = period2whole_vector(y, n);
    
    % y'*K*y + alpha*x'*x
    res = y' * semi_vector_multiplication(Ut, Vt, y);
    res = res + alpha * (x'*x);
    
    
end

