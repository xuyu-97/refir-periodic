function [return_filename] = data_generation(N, n, Lorder, poleub, noisetype, T, repi)
warning('off');
M = 1000 + N;

% phase = 2*pi*rand;
% omg = 2*pi/T;
% u = rand*sin(omg*(1:M)+phase)';
NumChannel = 1;
Period = T;
NumPeriod = M/Period;
u = idinput([Period,NumChannel,NumPeriod],'rgs');


t = (0:M-1)';

% [system, poles] = generate_linear_system_randomly(Lorder, poleub, 0.1);

% 
x1 = 0.9542;
x2 = 0.9758;
numerator = [-0.3377 0];
denominator = [1 -(x1+x2) x1*x2];
poles = [x1 x2];
system = tf(numerator,denominator,1);
% 
ytrue = lsim(system, u, t); 
[g,~] = impulse(system, t);


snr = 1;%10*(0.1+0.9*rand);
switch noisetype
    case 'gauss'
        noise = randn(N,1);
    case 'uniform'
        noise = rand(N,1)-0.5;
end
noise = (noise - mean(noise))*std(ytrue)/std(noise)/sqrt(snr);
y = ytrue(end-N+1:end) + noise;
u = u(end-N+1:end);
data = [u,y];

datainfo.data = data;
datainfo.ytrue = ytrue(end-N+1:end);
datainfo.LinearSystem = system;
datainfo.LinearSystemPoles = poles;
datainfo.LinearSystemImpulseResponse = g;
datainfo.snr = snr;
datainfo.noise = noise;
datainfo.LinearSysOrder = Lorder;
datainfo.NoiseVariance = var(noise);
datainfo.InputPeriod = T;


figure(1);
h = plot(0:n-1,g(1:n)); grid on;
title(['repi=' int2str(repi)]);
saveas(h,['Databank/Figs/system' int2str(repi) '.eps']);

return_filename = ['Databank/data_N' int2str(N) '_repi=' int2str(repi) '.mat'];
save(return_filename, 'datainfo');
end


function [system, poles] = generate_linear_system_randomly(order, poleub, polelb, fs_index)
% Randomly generate a stable lienar system
% the argument 'order' is the needed order of the system
if nargin < 4 || strcmp(fs_index,'fast')
    pole_max = 2; 
    pole_min = 0;
    value = 1;


   while pole_max > poleub || pole_min < polelb  || value > 1e6
       bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = 1; %bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 0;
        system = idpoly(md);
        
        poles = pole(system);
        pole_max = max(abs(poles));
        pole_min = min(abs(poles));
                
        value = max(abs([system.f system.b]));
   end

elseif strcmp(fs_index,'slow')

    pole_max = 0.5; % pole check
    value = 1;
    while pole_max < poleub || pole_max > 1 - eps || value > 1e6
        bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = 1;%bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 0;
        system = idpoly(md);
        pole_max = max(abs(pole(system)));
        value = max(abs([system.f system.b]));
    end
    poles = pole(system);
end

end
