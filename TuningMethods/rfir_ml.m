function [EstInfo] = rfir_ml(data, n, kernel)

u = data(:,1);
N = length(u);
Nd = N - n;
y = data(n+1:N,2);
Psi = CalculatePsi(u, n);
%QR decomposition of the unvaring part
Rd = triu(qr([Psi y]));
Rd = Rd(1:n+1,:);

%Hyper-parameter optimization
[lb,ub] = lubounds(kernel);

tic
hpini = ini_rfir(Rd, kernel, Nd);
time_ini = toc;

ff = @(x)nglglklhd(x, Rd, kernel, Nd);
options = optimoptions('fmincon', 'Display', 'off');
hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);
[cost, theta] = nglglklhd(hp, Rd, kernel, Nd);

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

function [hpini] = ini_rfir(Rd, kernel, Nd)
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
                    obj(nc,nlam,nrho,nsigsqr) = nglglklhd(hpstart, Rd, kernel, Nd);
                end
            end
        end
    end
    [~, indx] = min(obj(:));
    [indc, indlam, indrho, indsigsqr] = ind2sub([nc nlam nrho nsigsqr],indx);
    hpini = [c(indc) lam(indlam) rho(indrho) sigsqr(indsigsqr)]';
end

if strcmp(kernel, 'TC') || strcmp(kernel, 'SS')
    % hp = [c lam sigsqr]
    c = [-5 : 1 : 3]'; Lc = length(c);
    lam = [0.1:0.1:0.9]'; Llam = length(lam);
    sigsqr = [0.01: 0.03 :0.2]; Lsigsqr = length(sigsqr);

    obj = zeros(Lc,Llam, Lsigsqr);
    for nc = 1:Lc
        for nlam  = 1:Llam
                for nsigsqr = 1:Lsigsqr
                    hpstart = [c(nc) lam(nlam) sigsqr(nsigsqr)]';
                    obj(nc,nlam,nsigsqr) = nglglklhd(hpstart, Rd, kernel, Nd);
                end
        end
    end
    [~, indx] = min(obj(:));
    [indc, indlam, indsigsqr] = ind2sub([nc nlam nsigsqr],indx);
    hpini = [c(indc) lam(indlam) sigsqr(indsigsqr)]';
end

end


function [obj, theta] = nglglklhd(hyper, Rd, kernel, Nd)
warning('off', 'MATLAB:nearlySingularMatrix');

[~, nrd] = size(Rd);
n = nrd - 1;
sigsqr = hyper(end);

%Kernel selection

if strcmp(kernel, 'DC')
    % hp = [c lam rho sigsqr]
    c = exp(hyper(1));
    lam = hyper(2);
    rho = hyper(3);
    U = (lam*rho).^(1:n)';
    V = (lam/rho).^(1:n)';
    K = c*(tril(U*V')+triu(V*U',1));
end

if strcmp(kernel, 'TC')
    % hp = [c lam sigsqr]
    c = exp(hyper(1)); lam = hyper(2);
    %     U = c*(lam).^(1:n)'; V = ones(n,1);
    %     K = tril(U*V')+triu(V*U',1);
    t=lam.^(1:n);
    triuP = triu(repmat(t,n,1),1);
    K = c*(triuP + triuP' + diag(t));
end

if strcmp(kernel, 'SS')
    % hp = [c lam sigsqr]
    c = exp(hyper(1)); lam = hyper(2);
    u1 = -(lam^3).^(1:n)'/6; v1 = ones(n,1);
    u2 = (lam^2).^(1:n)'/2;  v2 = lam.^(1:n)';
    U = c*[u1 u2]; V = [v1 v2];
    K = tril(U*V')+triu(V*U',1);
end

chol_flag = 1;

% condition number check
% cond_flag = 1;

% if cond(K) > 1e+200
%     cond_flag = 0;
% end

try
     L = chol(K)';
catch
    try
        L = chol(K+eps*eye(n))';
    catch
        try
            chol_flag = 0;
            L = chol(K+1e-4*eye(n))';
        catch 
            obj = 1/eps;
            if nargout > 1
                theta = nan*ones(n,1);
            end
            return
        end
        
    end

end

if chol_flag == 0 %|| cond_flag == 0
    obj = 1/eps;
    if nargout > 1
        theta = nan*ones(n,1);
    end
else
    %QR decomposition of the varying part
    Rd1 = Rd(:,1:n); Rd2 = Rd(:,end);
    R = real(triu(qr([Rd1*L, Rd2; sqrt(sigsqr)*eye(n), zeros(n,1)])));
    R = R(1:n+1,:);
    R1 = R(1:n,1:n); R2 = R(1:n,end); r = R(end,end);
    clear R;

    %Empirical Bayes
    obj = (Nd-n)*log(sigsqr)+2*sum(log(abs(diag(R1))))+r^2/sigsqr;

    if nargout > 1
        R1inv = eye(n)/R1;
        theta = L*R1inv*R2;
    end
end
end


