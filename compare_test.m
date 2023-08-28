clear;

N = 1200; n = 500;  Nd = N - n; T = 5;
d = load(['Databank/data_N2000_repi=1.mat']);
data = d.datainfo.data(1:N,:);
u = data(:,1);
y = data(n+1:N,2);



hyper = [20;0.9000000;0.460000000000000]
% hyper = [-3;0.500000000000000;0.010000000000000];
c = exp(hyper(1)); lam = hyper(2);  sigsqr = hyper(3);


% %SS generators
% u1 = -(lam^3).^(1:n)'/6; v1 = ones(n,1);
% u2 = (lam^2).^(1:n)'/2;  v2 = lam.^(1:n)';
% U = c*[u1 u2]; V = [v1 v2];


%DC generators
% U = sqrt(c)*(lam*rho).^(1:n)'; V = sqrt(c)*(lam/rho).^(1:n)';

% TC generators
U = sqrt(c)*(lam).^(1:n)'; V = sqrt(c)*(1).^(1:n)';
% 
Psi = CalculatePsi(u, n);

K = tril(U*V')+triu(V*U',1);

% 
%  L = chol(K)';
%direct computation
quadr_part0 = y'*((Psi*K*Psi'+sigsqr*eye(Nd))\y);
logdet_part0 = logdet(Psi*K*Psi'+sigsqr*eye(Nd));
obj0 = quadr_part0 + logdet_part0;

tic
%regular case
Rd = triu(qr([Psi y]));
Rd = Rd(1:n+1,:);
% L = chol(K)'
try
    L = chol(K)';
catch
    L = chol(K+eps*eye(n))';
end
Rd1 = Rd(:,1:n); Rd2 = Rd(:,end);
R = real(triu(qr([Rd1*L, Rd2; sqrt(sigsqr)*eye(n), zeros(n,1)])));
R = R(1:n+1,:);
R1 = R(1:n,1:n); R2 = R(1:n,end); r = R(end,end);
quadr_part1 = r^2/sigsqr;
logdet_part1 = (Nd-n)*log(sigsqr)+2*sum(log(abs(diag(R1))));
obj1 = quadr_part1 + logdet_part1;
toc

%route1
tic

Psi_p = Psi(1:T,1:T);

quadr_part2 = method1_quadratic_form_inv(Psi_p, Nd, n, y, U', V', sigsqr);
logdet_part2 = method1_logdeterminant(Psi_p, Nd, n, U', V', sigsqr);
obj2 = quadr_part2 + logdet_part2;
toc

%route2
tic
Psi_p = Psi(1:T,1:T);

quadr_part3 = quadratic_form_inv(Psi_p, Nd, n, y, U', V', sigsqr);
logdet_part3 = logdeterminant(Psi_p, Nd, n, U', V', sigsqr);
obj3 = quadr_part3 + logdet_part3;
toc

err_q1 = abs(quadr_part0-quadr_part1)
err_q2 = abs(quadr_part0-quadr_part2)
err_q3 = abs(quadr_part0-quadr_part3)
err_l1 = abs(logdet_part0-logdet_part1)
err_l2 = abs(logdet_part0-logdet_part2)
err_l3 = abs(logdet_part0-logdet_part3)
