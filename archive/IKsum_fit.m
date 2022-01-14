clc
close all
clear variables


% raw trace & normalization with capacitance
load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
% load('ds_Ktrace_wt.mat')
% ds_Ktrace = ds_Ktrace_wt;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};

K_data = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx');
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

% X0 = [30.0, 7.0, 33.5, 7.0, ...
%       22.5, 7.7, 45.2, 2.058, 1200.0, ...
%       5.7, 2.058];
low_bd = [0.0, 2.0, 0.0, 2.0, ...
          0.0, 2.0, 0.0, 0.0, 170.0, ...
          2.0, 0.0];
upper_bd = [70.0, 50.0, 70.0, 50.0, ...
        35.0, 14.0, 80.0, 5.0, 5000.0, ...
        24.0, 5.0];
num_vars = 11;
fit_fn = @(X) Ktrace_fitness(X, ds_Ktrace, Iss);
options = optimoptions('ga','PlotFcn', @gaplotbestf, 'FitnessLimit',1000);
[param,fval] = ga(fit_fn,num_vars,[],[],[],[],low_bd,upper_bd,[],options);

