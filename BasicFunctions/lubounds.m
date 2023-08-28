function [lb,ub] = lubounds(kernel)
tol_lam = 0.01;
tol_sigsqr = eps;

if strcmp(kernel, 'DC')
% hp = [c lam rho sigsqr]
lb = [-inf; tol_lam; 0;  tol_sigsqr];
ub = [ inf; 1-tol_lam; 1; inf];
end

if strcmp(kernel, 'TC')  || strcmp(kernel, 'SS')
% hp = [c lam rho sigsqr]
lb = [-10; tol_lam;   tol_sigsqr];
ub = [ 10; 1-tol_lam; inf];
end
end
