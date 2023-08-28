function det_res = method1_logdeterminant(Phi_p, N, n, Ut, Vt, alpha)
% Calculate logdetH = logdet(Phi*K*Phi'+alpha*I_N)
    p = size(Phi_p, 1);
    % Phi = QR, obtain R: p*n
    [~, Rp] = HouseholderQR_periodic(Phi_p, N);
    
    % Compute R*K*R'+alpha*I_N (1:p,1:p)
    Rp = period2whole_matrix(Rp, p, n);
    Rp = Rp * (semi_matrix_multiplicaiton(Ut, Vt, Rp'));
    Rp = Rp + alpha * eye(p);
    
    % result
    det_res = (N-p)*log(alpha) + logdet(Rp);

end

