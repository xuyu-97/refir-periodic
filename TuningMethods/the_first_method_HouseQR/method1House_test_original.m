clear;
load('exm_Phi.mat');
Phi = Phi2;
[N, n] = size(Phi);
p = 70;
Phi_p = Phi(1:p, 1:p);
x = randn(N,1);
alpha = 0.9;

% Create instance of K=S(U,V) with rank k
k = 50;
%Ut = randn(k,n); Vt = randn(k,n);
% Ut = (0.9*0.8).^(1:n); Vt = (0.9/0.8).^(1:n);

c = exp(10); lam = 0.5;
u1 = -(lam^3).^(1:n)'/6; v1 = ones(n,1);
u2 = (lam^2).^(1:n)'/2;  v2 = lam.^(1:n)';
Ut = c*[u1 u2]'; Vt = [v1 v2]';
K = tril(Ut'*Vt) + triu(Vt'*Ut,1);

%% Comparison: x'*H*x 
% Reference
tic
H = Phi * K * Phi' + alpha * eye(N);
yref = x' * H * x;
toc

% Using method 2
tic
y = method1_quadratic_form(Phi_p, N, n, x, Ut, Vt, alpha);
toc

% Compare the accuracy
norm(yref - y)


%% Comparison: x'*inv(H)*x
% Reference
tic
yref = x'*inv(H)*x;
toc

% Using method 2
tic
y = method1_quadratic_form_inv(Phi_p, N, n, x, Ut, Vt, alpha);
toc

% Compare the accuracy
norm(yref - y)


%% Comparison: logdet(H)
% Reference
tic
det_ref = logdet(H);
toc

% Using method 2
tic
det_med = method1_logdeterminant(Phi_p, N, n, Ut, Vt, alpha);
toc

% Compare the accuracy
norm(det_ref - det_med)
