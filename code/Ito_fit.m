clc
close all
clear variables


%% Ito - KO
% import data
K_ko = readtable("./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx");
Ito_ko = K_ko.A3FF;

% voltage clamp protocol
holding_p = -70; %mV
holding_t = 4.5; %msec
P1 = 50; %mV
P1_t = 29.5; % msec
P2 = 50; % mV
P2_t = 29.5; % msec

% GA
num_vars = 6;
lower_bd = [0,0,0,0,0,0];
upper_bd = [60,60,60,60,60,60];
fit_fn = @(X) Ito_fitness(X,Ito_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);
[X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);

% Compare with the Rasmusson
[t,S,A,C] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
[t2,S2,A2,C2] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,[30,30,13.5,33.5,33.5,33.5]);

plot(t,A(:,61),'LineWidth',2)
hold on
plot(t2,A2(:,61),'LineWidth',2)
hold off
title('I_{to} Rapidly inactivating transient outward current - KO')
legend('GA','Rasmusson')
xlabel('sec')
ylabel('pA')

% Compare with the real data
fitted_Ito = A(:,61);
fitted_Ito(end)


%% Ito - WT
K_wt = readtable("./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx");
Ito_wt = K_wt.A3;

% voltage clamp protocol
holding_p = -70; %mV
holding_t = 4.5; %msec
P1 = 50; %mV
P1_t = 29.5; % msec
P2 = 50; % mV
P2_t = 29.5; % msec

% GA
num_vars = 6;
lower_bd = [0,0,0,0,0,0];
upper_bd = [60,60,60,60,60,60];
fit_fn = @(X) Ito_fitness(X,Ito_wt,holding_p,holding_t,P1,P1_t,P2,P2_t);
[X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);

% Compare with the Rasmusson
[t,S,A,C] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
[t2,S2,A2,C2] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,[30,30,13.5,33.5,33.5,33.5]);

plot(t,A(:,61),'LineWidth',2)
hold on
plot(t2,A2(:,61),'LineWidth',2)
hold off
title('I_{to} Rapidly inactivating transient outward current - WT')
legend('GA','Rasmusson')
xlabel('sec')
ylabel('pA')

% Compare with the real data
fitted_Ito = A(:,61);
fitted_Ito(end)
