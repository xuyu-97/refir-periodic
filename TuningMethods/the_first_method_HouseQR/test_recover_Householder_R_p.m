function R = test_recover_Householder_R_p(R_p, N, n)
% R_p is obtained by HouseholderQR_periodic(Phi_p, N).
% Phi_p is the periodic part of N*n Phi.
% Recover R_p to N*n R
    p = size(R_p, 1);
    R = zeros(N, n);
    R(1:p, :) = period2whole_matrix(R_p, p, n);
end

