function res = quadratic_form_inv(Phi_p, N, n, x, Ut, Vt, alpha)
% Input: period part Phi_p, dimension of Phi: N*n, period p,
%        N-dim x, semi-separable representation U & V, and alpha
% H = Phi * K * Phi' + alpha * I
% Phi is periodic Toeplitz matrix with period p. K = S(U,V) with rank k.
% Output: x'* inv(H) * x
    p = size(Phi_p, 1);
    
    % Phi = QR, obtain R: p*n
    R = period_Toeplitz_qr_getR(Phi_p, N);
    R = period2whole_matrix(R, p, n);
    
    % 1) (fixed) y=Phi'*x
    y = Phi_p_times_x(Phi_p', N, x);
    y = period2whole_vector(y, n);
    
    % 2) (fixed) x'*x
    le = x'*x;
    
    % y'*Hbar*y
    % 3) K'*y = K*y and RK'y
    % 4) y'Ky
    z1 = semi_vector_multiplication(Ut, Vt, y);
    z2 = y' * z1;
    z1 = R * z1;
    
    % 5) z1'*(R*K*R'+alpha*I_p)^(-1) * z1
    Inside = semi_matrix_multiplicaiton(Ut,Vt,R');
    Inside = R * Inside;
    Inside = Inside + alpha*eye(p);
    z1 = z1' * (Inside \ z1);
    
    % 6) (z2 - z1)/alpha
    z = (z2-z1)/alpha;
    
    % Result
    res = (le - z)/alpha;
    
end

