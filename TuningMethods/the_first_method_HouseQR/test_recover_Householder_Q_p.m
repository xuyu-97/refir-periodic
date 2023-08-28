function Q = test_recover_Householder_Q_p(Q_p, N)
% Input: (p+2)*p Q_p obtained in HouseholderQR_periodic.m and 
% Phi with size N*p
% Output: recover Q: N*N
    p = size(Q_p, 2);
    Q = eye(N);
    
    % P_k = [I_{k-1}, 0; 0, I_{N-(k-1)} - b_k * v_k * v_k']
    % Q = P_1 * P_2 * ... * P_p
    for k = 1:p
        b_k = Q_p(p+2,k);
        v_k_p = Q_p(1:(p+1),k);
        v_k = [v_k_p(1); period2whole_vector(v_k_p(2:(p+1)), N-k)];
        P_k = [eye(k-1), zeros(k-1,N-k+1); ...
               zeros(N-k+1,k-1), eye(N-k+1)-b_k*v_k*v_k'];
        Q = Q * P_k;
    end
end

