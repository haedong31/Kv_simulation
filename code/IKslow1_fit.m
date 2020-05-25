clc
close all
clear variables


K_ko = readtable("./potassium-KO.xlsx");
K_wt = readtable("./potassium-WT.xlsx");
IKslow1_ko = K_ko.A2FF;
IKslow1_wt = K_wt.A2;
tau1_ko = K_ko.Tau2FF;
tau1_wt = K_wt.Tau2;
Y_ko = table(IKslow1_ko, tau1_ko);
Y_wt = table(IKslow1_wt, tau1_wt);
Y_ko.Properties.VariableNames = {'AMP', 'TAU'};
Y_wt.Properties.VariableNames = {'AMP', 'TAU'};
% X = [22.5, 7.7, 45.2, 5.7, 45.2, 5.7];

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 450; % ms
P1 = 50; % mV
P1_t = 25 * 1000; % ms
P2 = 50; % mV
P2_t = P1_t; % ms

% run GA
num_vars = 6;
fit_fn_ko = @(X) IKslow1_fitness(X,Y_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);
[params_ko, fval] = ga(fit_fn_ko,num_vars);

fit_fn_wt = @(X) IKslow1_fitness(X,Y_wt,IKslow1_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);
[params_wt, fval] = ga(fit_fn_wt,num_vars);
