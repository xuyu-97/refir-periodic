function [Q_p, R_p] = HouseholderQR_periodic(Phi_p, N)
% Phi is a N*p periodic matrix with period p
% Input: periodic part Phi_p, and row size N
% Output: periodic part of Householder vectors Q; R
    
    p = size(Phi_p, 1);
    Q_p = zeros(p+2, p); % Q_p(1:p+1,:) for v_p, Q_p(p+2,:) for b
    R_p = zeros(p, p);
    
    for k = 1:p
    % Calculate Householder vectors and b to obtain Q_p(:, k)
    % P_k = [I_{k-1}, 0; 0, I_{N-(k-1)} - b_k * v_k * v_k']
        [v_p, b] = house_period(Phi_p(:,1), N-(k-1));
        Q_p(1:p+1, k) = v_p;
        Q_p(p+2, k) = b;
    
    % Update Phi_p (Phi)
%         Phi_p_1 = [Phi_p(2:p,:); Phi_p(1,:)];
%         m = floor((N-k)/p); r = mod(N-k, p);
%         Vp = v_p(1)*Phi_p(1,:) + m*v_p(2:(p+1))'*Phi_p_1 + v_p(2:(r+1))'*Phi_p_1(1:r,:);
%         Vp = b * v_p * Vp;
%         Phi_p = [Phi_p(1,:)-Vp(1,:); Phi_p_1 - Vp(2:(p+1),:)];
%         R_p(k,k:p) = Phi_p(1,:);
%         Phi_p = Phi_p(2:(p+1), 2:(p-k+1));
        
        
        Phi_p = [Phi_p(2:p,:); Phi_p(1,:)];
        m = floor((N-k)/p); r = mod(N-k, p);
        Vp = v_p(1)*Phi_p(p,:) + m*v_p(2:(p+1))'*Phi_p + v_p(2:(r+1))'*Phi_p(1:r,:);
        Vp = b * v_p * Vp;
        Phi_p = [Phi_p(p,:)-Vp(1,:); Phi_p - Vp(2:(p+1),:)];
        R_p(k,k:p) = Phi_p(1,:);
        Phi_p = Phi_p(2:(p+1), 2:(p-k+1));
    end
    
end



function [v_p, b] = house_period(x_p, n)
    % Given a size-n vector x with period p, calculate householder vector's
    % period part v_p and b, such that P = I_n - bvv'
    
    p = size(x_p, 1);
    % sigma = x(2:n)'*x(2:n)
    x_p_1 = [x_p(2:p); x_p(1)];
    m = floor((n-1)/p);  r = mod(n-1, p);
    sigma = m * x_p_1' * x_p_1 + x_p_1(1:r)'*x_p_1(1:r);
    
    % v = [1; x(2:n)], just consider v_p: v(1:p+1) part
    v_p = [1; x_p_1];
    
    if sigma == 0
        b = 0;
    else
        mu = sqrt(x_p(1)^2 + sigma);
        if x_p(1) <= 0
            v_p(1) = x_p(1) - mu;
        else
            v_p(1) = -sigma / (x_p(1) + mu);
        end
        b = (2*v_p(1)^2) / (sigma + v_p(1)^2);
        v_p = v_p / v_p(1);
    end
end