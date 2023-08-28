function [Rp, x] = method1simu_fix_part(Phi_p, N, n, x)
% Input: period part Phi_p, dimension of Phi: N*n. Phi has period p.
% Purpose: H = Phi * K * Phi' + alpha * I
% We want to compute x'*inv(H)*x + logdet(H)
% We use the first route: Phi = QR (full)
% -- x'*inv(H)*x = x'*Q*(R*K*R'+alpha*I)*Q'*x
% -- logdet(H)   = (N-p)*log(alpha) + logdet(R*K*R'+alpha*I)
% Output: fixed part R and Q'*x

    p = size(Phi_p, 1);
    [Qp, Rp] = HouseholderQR_periodic(Phi_p, N);
    
    % Q' * x = Pp * ... * P1 * x
    for k = 1:p
        b_k = Qp(p+2,k);
        v_k = Qp(1:(p+1),k);
        v_k = [v_k(1); period2whole_vector(v_k(2:(p+1)), N-k)];
        x(k:N) = x(k:N) - b_k*v_k*(v_k'*x(k:N));
    end
    
    % extend Rp into a periodic p*n matrix
    Rp = period2whole_matrix(Rp, p, n);
end

