function res = method1_quadratic_form_inv(Phi_p, N, n, x, Ut, Vt, alpha)
%  Input: period part Phi_p, dimension of Phi: N*n, period p,
%        N-dim x, semi-separable representation U & V, and alpha
% H = Phi * K * Phi' + alpha * I
% Phi is periodic Toeplitz matrix with period p. K = S(U,V) with rank k.
% Output: x'* inv(H) * x
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
    
    % Compute R*K*R'+alpha*I_N (1:p,1:p)
    Rp = period2whole_matrix(Rp, p, n);
    Rp = Rp * (semi_matrix_multiplicaiton(Ut, Vt, Rp'));
    Rp = Rp + alpha * eye(p);
    
    % Result
    res = x(1:p)' * (Rp\x(1:p));
    res = res + (x((p+1):N)' * x((p+1):N))/alpha;
end

