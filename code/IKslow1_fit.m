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
t = 0:1:25*1000;
V = zeros(1, length(t));
V(1:450) = -70;
V(452:end) = 50;

options = optimoptions('fminunc', 'Algorithm','quasi-newton', 'Display','iter');
problem.options = options;
problem.x0 = [33.5, 7.0, 22.5, 7.0, 45.2, 5.7, 1050.0, 45.2, 2.058];
problem.objective = @(X)IKslow1_fitness(X, Y_ko, t, V);
problem.solver = 'fminunc';
est_x = fminunc(problem);
