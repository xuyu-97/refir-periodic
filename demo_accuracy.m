clear;

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);
addpath([pwd '/TuningMethods/the_first_method_HouseQR']);
addpath([pwd '/TuningMethods/the_second_method']);

%Specification
Nmax = 2000;
N = 600;

Lorder = 10;
kernel = {'TC'};

%period of input signal
T = 40;
n = 50;
MaxRepi = 80;

%Data generation
snr = 10;
for repi = 1:MaxRepi
    data_generation(Nmax, n, Lorder, 0.9, 'gauss', T, repi); %pole <= 0.8 works for route2
end

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

for repi = 1:MaxRepi
    fprintf('-----------------repi = %i-----------------\n', repi);
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    ytrue = d.datainfo.ytrue(n+1:N);
    gtrue = d.datainfo.LinearSystemImpulseResponse(2:n+1);

    %Regular ReFir method
    EstInfo = rfir_ml(data, n, kernel);
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
%     EstInfo2 = rfir_ml_route2_new(data, n, kernel, T);
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

%Store the results
Result_regular.HP = HP;
Result_regular.HPini = HPi;
Result_regular.EFIT = EFIT;
Result_regular.GFIT = GFIT;
Result_regular.COST = COST;


Result_route1.HP = HP1;
Result_route1.HPini = HPi1;
Result_route1.EFIT = EFIT1;
Result_route1.GFIT = GFIT1;
Result_route1.COST = COST1;


% Result_route2.HP = HP2;
% Result_route2.HPini = HPi2;
% Result_route2.EFIT = EFIT2;
% Result_route2.GFIT = GFIT2;
% Result_route2.COST = COST2;

save('Result_regular.mat', 'Result_regular');
save('Result_route1.mat', 'Result_route1');
% save('Result_route2.mat', 'Result_route2');

% figure(1)
% p_boxplot([GFIT' GFIT1' ],0,100,{'RFIR','RFIR-r1'},'GFIT','')
% figure(2)
% p_boxplot([EFIT' EFIT1'],0,100,{'RFIR','RFIR-r1'},'EFIT','')

p_boxplot(GFIT'- GFIT1',-1e-5,1e-5,{''},'','')

format long
mean([GFIT' GFIT1'])

