clc
close all
clear variables


K_data = readtable('./potassium-KO.xlsx');
% K_data = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx');

Iss = K_data.IssFF;
% Iss = K_data.Iss;
Iss = nanmean(Iss);

Ito = K_data.A3FF;
% Ito = K_data.A3;
Ito = nanmean(Ito);

tau_to = K_data.Tau3FF;
% tau_to = K_data.Tau3;
tau_to = nanmean(tau_to);

IKslow1 = K_data.A2FF;
% IKslow1 = K_data.A2;
IKslow1 = nanmean(IKslow1);

tau1 = K_data.Tau2FF;
% tau1 = K_data.Tau2;
tau1 = nanmean(tau1);

IKslow2 = K_data.A1FF;
% IKslow2 = K_data.A1;
IKslow2 = nanmean(IKslow2);

tau2 = K_data.Tau1FF;
% tau2 = K_data.Tau1;
tau2 = nanmean(tau2);

cap = K_data.CapFF;
% cap = K_data.Cap;
cap = nanmean(cap);

X0 = [7.7, 45.2, 5.7, 1200.0, 45.2];
low_bd = [2.0, 10.0, 2.0, 200.0, 0.0];

num_vars = 5;
fit_fn = @(X) IKslow2_fitness(X, IKslow2, tau2);
opts = optimoptions('ga', 'PlotFcn','gaplotbestf', 'PopulationSize',40, 'MaxTime',18000);
[param,fval] = ga(fit_fn,num_vars,[],[],[],[],low_bd,[],[],opts);
