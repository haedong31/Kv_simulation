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

%% I_to - WT
% K_wt = readtable("./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx");
