function [y, R, le] = method2simu_fix(Phi_p, N, n, x)
% Input: period part Phi_p, dimension of Phi: N*n. Phi has period p.
% Purpose: H = Phi * K * Phi' + alpha * I
% We want to compute x'*inv(H)*x + logdet(H)
% We use the second route: Phi = QR (thin)
% Output: fixed part
    p = size(Phi_p, 1);
    
    % Phi = QR, obtain R: p*n
    R = period_Toeplitz_qr_getR(Phi_p, N);
    R = period2whole_matrix(R, p, n);
    
    % 1) (fixed) y=Phi'*x
    y = Phi_p_times_x(Phi_p', N, x);
    y = period2whole_vector(y, n);

    % 2) (fixed) x'*x
    le = x'*x;

end

