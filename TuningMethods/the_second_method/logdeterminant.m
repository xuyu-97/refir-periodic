function det_res = logdeterminant(Phi_p, N, n, Ut, Vt, alpha)
% Calculate det(Phi*K*Phi'+alpha*I_N)

    p = size(Phi_p, 1);
    % Phi = QR, obtain R: p*n
    R = period_Toeplitz_qr_getR(Phi_p, N);
    R = period2whole_matrix(R, p, n);
    
    % alpha*I_p+RKR'
%     Inside = semi_matrix_multiplicaiton(Ut,Vt,R');
%     Inside = R * Inside;
%     Inside = Inside + alpha*eye(p);
    R = R * (semi_matrix_multiplicaiton(Ut,Vt,R'));
    R = R + alpha*eye(p);
    
    % result
    det_res = (N-p)*log(alpha) + logdet(R);
end

