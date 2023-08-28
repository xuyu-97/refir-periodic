clear;
% % 
% c=parcluster('local');
% c.NumWorkers= 3;
% parpool(3);
addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);
addpath([pwd '/TuningMethods/the_first_method_HouseQR']);
addpath([pwd '/TuningMethods/the_second_method']);


%Specification
Nmax = 10000;
N = 10000;
nmax = 4800;
Lorder = 10;
kernel = {'TC'};

%period of input signal
T = 200;
nrange =   [300, 600, 1200 ,2400, 4800];
MaxRepi = length(nrange);

%Data generation
snr = 10;
% data_generation(Nmax, nmax, Lorder, 0.95, 'gauss', T, 1);

EFIT = [];
GFIT = [];
HP = [];
HPi = [];
COST = [];

EFIT1 = [];
GFIT1 = [];
HP1 = [];
HPi1 = [];
COST1 = [];

EFIT2 = [];
GFIT2 = [];
HP2 = [];
HPi2 = [];
COST2 = [];

t = zeros(MaxRepi,1);
t1 = zeros(MaxRepi,1);
t2 = zeros(MaxRepi,1);

avcount = 10;

for repi =1:MaxRepi
    n = nrange(repi);
    fprintf('-----------------n = %i-----------------\n',n);
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(1) '.mat']);
    data = d.datainfo.data(1:N,:);
    ytrue = d.datainfo.ytrue(n+1:N);
    gtrue = d.datainfo.LinearSystemImpulseResponse(2:n+1);
    
    tt = 0; tt1 = 0; tt2 = 0;
    for iter = 1:avcount
    fprintf('iter = %i\n',iter);
    
    %Regular ReFir method
    EstInfo = rfir_ml(data, n, kernel);
    tt = tt + EstInfo.time_ini;
    hp = EstInfo.hp;
    cost = EstInfo.cost;
    ghat = EstInfo.ghat;
    yhat = EstInfo.yhat;
    efit= gof(ytrue,yhat);
    gfit = gof(gtrue,ghat);
    EFIT = [EFIT efit];
    GFIT = [GFIT gfit];
    COST = [COST cost];
    HP = [HP hp];
    HPi = [HPi EstInfo.hpini];

    %ReFir method by route1
    EstInfo1 = rfir_ml_route1_new(data, n, kernel, T);
    tt1 = tt1 + EstInfo1.time_ini;
    hp1 = EstInfo1.hp;
    cost1 = EstInfo1.cost;
    ghat1 = EstInfo1.ghat;
    yhat1 = EstInfo1.yhat;
    efit1 = gof(ytrue,yhat1);
    gfit1 = gof(gtrue,ghat1);
    EFIT1 = [EFIT1 efit1];
    GFIT1 = [GFIT1 gfit1];
    COST1 = [COST1 cost1];
    HP1 = [HP1 hp1];
    HPi1 = [HPi1 EstInfo1.hpini];

%     %ReFir method by route2
%     
%     EstInfo2 = rfir_ml_route2_new(data, n, kernel, T);
%     tt2 = tt2 + EstInfo2.time_ini;
%     hp2 = EstInfo2.hp;
%     cost2 = EstInfo2.cost;
%     ghat2 = EstInfo2.ghat;
%     yhat2 = EstInfo2.yhat;
%     efit2 = gof(ytrue,yhat2);
%     gfit2 = gof(gtrue,ghat2);
%     EFIT2 = [EFIT2 efit2];
%     GFIT2 = [GFIT2 gfit2];
%     COST2 = [COST2 cost2];
%     HP2 = [HP2 hp2];
%     HPi2 = [HPi2 EstInfo2.hpini];
    end
    t(repi) = tt/avcount;
    t1(repi) = tt1/avcount;
%     t2(repi) = tt2/avcount;
end

%Store the results
Result_regular.HP = HP;
Result_regular.HPini = HPi;
Result_regular.EFIT = EFIT;
Result_regular.GFIT = GFIT;
Result_regular.COST = COST;
Result_regular.time = t;

Result_route1.HP = HP1;
Result_route1.HPini = HPi1;
Result_route1.EFIT = EFIT1;
Result_route1.GFIT = GFIT1;
Result_route1.COST = COST1;
Result_route1.time = t1;

% Result_route2.HP = HP2;
% Result_route2.HPini = HPi2;
% Result_route2.EFIT = EFIT2;
% Result_route2.GFIT = GFIT2;
% Result_route2.COST = COST2;
% Result_route2.time = t2;

save('Result_regular.mat', 'Result_regular');
save('Result_route1.mat', 'Result_route1');
% save('Result_route2.mat', 'Result_route2');

% % 3 routes logplot
% figure(1)
% semilogy(100:100:2200,t(1:22)); hold on; 
% semilogy(100:100:2200,t1(1:22),'-*'); hold on; 
% semilogy(100:100:2200,t2(1:22),'--'); hold on; grid on;
% legend('RFIR','RFIR-r1','RFIR-r2');
% xlabel('n'); ylabel('averaged computation time (s)');
% 
% % 3 routes plot
% figure(2)
% plot(100:100:2200,t(1:22)); hold on; 
% plot(100:100:2200,t1(1:22),'-*'); hold on; 
% plot(100:100:2200,t2(1:22),'--'); hold on; grid on;
% legend('RFIR','RFIR-r1','RFIR-r2');
% xlabel('n'); ylabel('averaged computation time (s)');

% 
% % 2 routes logplot
% figure(3)
% semilogy(100:100:2200,t(1:22)); hold on; 
% semilogy(100:100:2200,t1(1:22),'-*'); hold on; grid on;
% legend('RFIR','RFIR-hss');
% xlabel('n'); ylabel('averaged computation time (s)');
% 
% % 2 routes logplot
% figure(4)
% plot(100:100:2200,t(1:22)); hold on; 
% plot(100:100:2200,t1(1:22),'-*'); hold on; grid on;
% legend('RFIR','RFIR-hss');
% xlabel('n'); ylabel('averaged computation time (s)');






