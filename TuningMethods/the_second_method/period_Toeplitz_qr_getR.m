function R = period_Toeplitz_qr_getR(Phi_p, N)
% Input: period part of period Toeplitz matrix Phi(1:p,1:p) with size N*p
% Output: R where Phi=QR
% Algorithm: Qiao, 1988

    p = size(Phi_p, 2);
    r = mod(N, p);
    x_p = [Phi_p(2:p,1); Phi_p(1,1)]; % size of x: N-1. x_p=x(1:p)
    x_R = flip([Phi_p(2+r:p,1); Phi_p(1:r,1)]); % size p-1
    y = Phi_p(1, 2:p); 
    
% W: Givens transformation from 1 to n-1
% W_i: n-dim, on (i, n) rows
% The first row is cos (cw_i), the second row is sin (sw_i).
    W_rot = zeros(2, p-1);    
    
% U, V: Givens transformation from 1 to n-1, n-dim, on (i,i+1) rows
    U_rot = zeros(2, p-1);
    V_rot = zeros(2, p-1);
                                
% Initialize R, R(1,1) and eta
    R = zeros(p);
    temp = Phi_p(:,1).^2;
    R(1,1) = floor(N/p)*sum(temp) + sum(temp(1:r));
    R(1,1) = sqrt(R(1,1));
    eta = R(1,1);

% R_1: T_1 = Q_1 * R_1, both are p-1 by p-1
    R_1 = zeros(p-1);
% R_bar and a for Stage 2 K_2 = [Phi(1,1), y'; a, R_1]
    R_bar = zeros(p-1);        % (p-1)dim matrix
    a = zeros(p-1, 1);
    K2_1st_col = [Phi_p(1,1); a];
    
% Calculate T_1' * x in (2.8) such that R_1' * a = Phi_1' * x using FFT
    Phi_1_t_times_x = Phi_p_times_x_p(Phi_p', N-1, x_p);
    Phi_1_t_times_x = Phi_1_t_times_x(1:(p-1));
    
    
% ======================================================================= %
% QR decomposition: Get R, Stage 1 and Stage 2
    for k=1:(p-1)
    % Stage 1: {R(:,k), W_1, ..., W_{k-1} -> R_1(:,k), W_k}
        t = x_R(k);
        for i = 1:(k-1)
            R_1(i,k) = (R(i,k) - W_rot(2,i) * t) / W_rot(1,i);
            t = -W_rot(2,i) * R_1(i,k) + W_rot(1,i) * t;
        end
        R_1(k,k) = sqrt(R(k,k)^2 - t^2);
        W_rot(1,k) = R_1(k,k) / R(k,k);
        W_rot(2,k) = t / R(k,k);
        
    % Stage 2: {R_1(:,k), V_1, ..., V_{k-1}, U_1, ..., U_{k-1}}
    % -> R(;,k+1), U_k, V_k
       
    % Determine R_bar(:,k) and V_k
        R_bar(1,k) = y(k);
        for i=1:(k-1)
            [cvi, svi] = deal(V_rot(1,i), V_rot(2,i));
            [R_bar(i,k), R_bar(i+1,k)] = deal(cvi*R_bar(i,k) + svi*R_1(i,k),...
                                             -svi*R_bar(i,k) + cvi*R_1(i,k));
        end
        ga = sqrt(R_bar(k,k)^2 + R_1(k,k)^2);
        V_rot(1,k) = R_bar(k,k) / ga;
        V_rot(2,k) = R_1(k,k) / ga;
        R_bar(k,k) = ga;
        
    % Solve a(k), we have computed Phi_1'*x
        a(k) = (Phi_1_t_times_x(k) - R_1(1:(k-1),k)' * a(1:(k-1))) / R_1(k,k);
        K2_1st_col(k+1) = a(k);
        
    % Premultiply V_k to V_{k-1}...V_1 * K2_1st_col
        [K2_1st_col(k), K2_1st_col(k+1)] = ...
            deal(V_rot(1,k) * K2_1st_col(k) + V_rot(2,k) * K2_1st_col(k+1),...
                 -V_rot(2,k) * K2_1st_col(k) + V_rot(1,k) * K2_1st_col(k+1));
        
    % Compute U_k and eta_{k+1}
        eta_n = sqrt(eta^2 - K2_1st_col(k)^2);
        U_rot(1,k) = K2_1st_col(k) / eta;
        U_rot(2,k) = eta_n / eta;
        eta = eta_n;


    % Get R(:,k+1)
        for i = k:-1:1
            if i == p-1
                [R_bar(i,k), R(i+1,k+1)] = ...
                deal(U_rot(1,i) * R_bar(i,k), -U_rot(2,i) * R_bar(i,k));
            else
                [R_bar(i,k), R(i+1,k+1)] = ...
                deal(U_rot(1,i) * R_bar(i,k) + U_rot(2,i) * R_bar(i+1,k),...
                     -U_rot(2,i) * R_bar(i,k) + U_rot(1,i) * R_bar(i+1,k));
            end   
        end
        
        R(1,k+1) = R_bar(1,k);
    end

end

