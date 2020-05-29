clc
close all
clear variables


% raw trace & normalization with capacitance
load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
% load('ds_Ktrace_wt.mat')
% ds_Ktrace = ds_Ktrace_wt;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};

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

ds_Ktrace.I = ds_Ktrace.I ./ cap;

num_vars = 22;
fit_fn = @(X) Ktrace_fitness(X, ds_Ktrace, Iss, Ito, IKslow1, IKslow2, tau_to, tau1, tau2);
options = optimoptions('ga','PlotFcn', @gaplotbestf, 'FitnessLimit',1500);
[param,fval] = ga(fit_fn,num_vars,options);

plot(ds_Ktrace.time, ds_Ktrace.I, 'LineWidth',2);
hold on
plot(t, IKsum, 'LineWidth',2);
hold off
legend('Raw Trace', 'Simulated')
