function res = method1simu_obj(N, R, Qtx, Ut, Vt, alpha)
% Purpose: H = Phi * K * Phi' + alpha * I
% We want to compute x'*inv(H)*x + logdet(H)
% We use the first route: Phi = QR (full), 
% -- x'*inv(H)*x = x'*Q*(R*K*R'+alpha*I)*Q'*x
% -- logdet(H)   = (N-p)*log(alpha) + logdet(R*K*R'+alpha*I)

% Input: fixed part R, Qtx(Q'*x), N; change part Ut, Vt, alpha 
% Output: objective function x'*inv(H)*x + logdet(H)

    p = size(R,1);
    
    % Compute R*K*R'+alpha*I_N (1:p,1:p)
    
    R = R * (semi_matrix_multiplicaiton(Ut, Vt, R'));
    R = R + alpha * eye(p);
    
%     K = tril(Ut'*Vt) + triu(Vt'*Ut,1);
%     R = R*K*R' + alpha * eye(p);
    
    % result: x'*inv(H)*x + logdet(H)
    res_inv = Qtx(1:p)' * (R\Qtx(1:p));
    res_inv = res_inv + (Qtx((p+1):N)' * Qtx((p+1):N))/alpha;
    res_logdet = (N-p)*log(alpha) + logdet(R);
    res = res_inv + res_logdet; 
end

