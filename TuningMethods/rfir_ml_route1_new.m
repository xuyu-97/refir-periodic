function [EstInfo] = rfir_ml_route1_new(data, n, kernel, T)

u = data(:,1);
N = length(u);

y = data(n+1:N,2);
Nd = length(y);
Psi = CalculatePsi(u, n);
%Periodic part of Psi
Psi_p = Psi(1:T,1:T);

%Hyper-parameter optimization
[lb,ub] = lubounds(kernel);
[R, Qty] = method1simu_fix_part(Psi_p, Nd, n, y);

tic
hpini = ini_rfir(Qty, R, kernel, n);
time_ini = toc;

ff = @(x)nglglklhd(x, Qty, R, kernel, n);
options = optimoptions('fmincon', 'Display', 'off');
hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);
cost = nglglklhd(hp, Qty, R, kernel, n);

%Calculate the estimates with the optimal hyper-parameters
if strcmp(kernel, 'DC')
    c = exp(hp(1)); lam = hp(2); rho = hp(3);
    U = sqrt(c)*(lam*rho).^(1:n)'; V = sqrt(c)*(lam/rho).^(1:n)';
    K = tril(U*V')+triu(V*U',1);
end


if strcmp(kernel, 'TC')
    c = exp(hp(1)); lam = hp(2);
    %     U = c*lam.^(1:n)'; V = ones(n,1);
    %     K = tril(U*V')+triu(V*U',1);
    t=lam.^(1:n);
    triuP = triu(repmat(t,n,1),1);
    K = c*(triuP + triuP' + diag(t));
end

if strcmp(kernel, 'SS')
    % hp = [c lam sigsqr]
    c = exp(hp(1)); lam = hp(2);
    u1 = -(lam^3).^(1:n)'/6; v1 = ones(n,1);
    u2 = (lam^2).^(1:n)'/2;  v2 = lam.^(1:n)';
    U = c*[u1 u2]; V = [v1 v2];
    K = tril(U*V')+triu(V*U',1);
end

sigsqr = hp(end);

try
    L = chol(K)';
    R = real(triu(qr([Psi*L, y; sqrt(sigsqr)*eye(n), zeros(n,1)])));
    R = R(1:n+1,:);
    R1 = R(1:n,1:n); R2 = R(1:n,end);
    R1inv = eye(n)/R1;
    theta = L*R1inv*R2;
catch
    theta = K*Psi'*((Psi*K*Psi'+sigsqr*eye(N-n))\y);
end


%Estimation information collection
yhat = Psi*theta;
EstInfo.ghat = theta;
EstInfo.hpini = hpini;
EstInfo.time_ini = time_ini;
EstInfo.hp = hp;
EstInfo.kernel = kernel;
EstInfo.yhat = yhat;
EstInfo.Psi = Psi;
EstInfo.cost = cost;
end

function [hpini] = ini_rfir(Qty, R, kernel, n)
%Optimization by grid search or other methods
if strcmp(kernel, 'DC')
    % hp = [c lam rho sigsqr]
    c = [-15 : 1 : 5]'; Lc = length(c);
    lam = [0.1 0.2 0.3:0.05:0.9]'; Llam = length(lam);
    rho = [0.25:0.15:0.95]'; Lrho = length(rho);
    sigsqr = [0.01: 1 :10]; Lsigsqr = length(sigsqr);

    obj = zeros(Lc,Llam,Lrho, Lsigsqr);
    for nc = 1:Lc
        for nlam  = 1:Llam
            for nrho = 1:Lrho
                for nsigsqr = 1:Lsigsqr
                    hpstart = [c(nc) lam(nlam) rho(nrho) sigsqr(nsigsqr)]';
                    obj(nc,nlam,nrho,nsigsqr) = nglglklhd(hpstart, Qty, R, kernel, n);
                end
            end
        end
    end
    [~, indx] = min(obj(:));
    [indc, indlam, indrho, indsigsqr] = ind2sub([nc nlam nrho nsigsqr],indx);
    hpini = [c(indc) lam(indlam) rho(indrho) sigsqr(indsigsqr)]';
end

if strcmp(kernel, 'TC')  || strcmp(kernel, 'SS')
    % hp = [c lam sigsqr]
    c = [-5 : 1 : 3]'; Lc = length(c);
    lam = [0.1:0.1:0.9]'; Llam = length(lam);
    sigsqr = [0.01: 0.03 :0.2]; Lsigsqr = length(sigsqr);


    obj = zeros(Lc,Llam, Lsigsqr);
    for nc = 1:Lc
        for nlam  = 1:Llam
                for nsigsqr = 1:Lsigsqr
                    hpstart = [c(nc) lam(nlam) sigsqr(nsigsqr)]';
                    obj(nc,nlam,nsigsqr) = nglglklhd(hpstart, Qty, R, kernel, n);
                end
        end
    end
    [~, indx] = min(obj(:));
    [indc, indlam, indsigsqr] = ind2sub([nc nlam nsigsqr],indx);
    hpini = [c(indc) lam(indlam) sigsqr(indsigsqr)]';
end

end


function [obj] = nglglklhd(hyper, Qty, R, kernel, n)
warning('off', 'MATLAB:nearlySingularMatrix');
Nd = length(Qty);

sigsqr = hyper(end);

%Kernel selection

if strcmp(kernel, 'DC')
    % hp = [c lam rho sigsqr]
    c = exp(hyper(1));
    lam = hyper(2);
    rho = hyper(3);
    U = sqrt(c)*(lam*rho).^(1:n)';
    V = sqrt(c)*(lam/rho).^(1:n)';
end

if strcmp(kernel, 'TC')
    c = exp(hyper(1)); lam = hyper(2);
    U = c*lam.^(1:n)'; V = ones(n,1);
end

if strcmp(kernel, 'SS')
    % hp = [c lam sigsqr]
    c = exp(hyper(1)); lam = hyper(2);
    u1 = -(lam^3).^(1:n)'/6; v1 = ones(n,1);
    u2 = (lam^2).^(1:n)'/2;  v2 = lam.^(1:n)';
    U = c*[u1 u2]; V = [v1 v2];
end
%Empirical Bayes
Ut = U'; Vt = V';
% obj = method1_quadratic_form_inv(Psi_p, Nd, n, y, Ut, Vt, sigsqr) ...
%     + method1_logdeterminant(Psi_p, Nd, n, Ut, Vt, sigsqr);
obj = method1simu_obj(Nd, R, Qty, Ut, Vt, sigsqr);
end


