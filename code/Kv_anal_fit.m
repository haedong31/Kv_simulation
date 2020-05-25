clc
close all
clear variables


load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};

K_data = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx');

Iss = K_data.IssFF;
Iss = nanmean(Iss);

cap = K_data.CapFF;
cap = nanmean(cap);

ds_Ktrace.I = ds_Ktrace.I ./ cap;

options = optimoptions('fminunc', 'Algorithm','quasi-newton', 'Display','iter');
problem.options = options;
problem.x0 = [-22.5, 7.7, 45.2, 5.7, 33.5, 7.0, ...
    33.5, 7.0, -22.5, 7.0, 45.2, 5.7, 1050.0, 45.2, 2.058,...
    33.5, 7.0, -22.5, 7.0, 45.2, 5.7, 1050.0, 45.2, 2.058];
problem.objective = @(X)Kv_anal_fitness(X, ds_Ktrace, Iss);
problem.solver = 'fminunc';
est_x = fminunc(problem);

[~, peak_idx] = max(ds_Ktrace.I);
Iksum = Kv_anal(ds_Ktrace.time(peak_idx:end), 50, est_x);
plot(ds_Ktrace.time, ds_Ktrace.I)
hold on
plot(ds_Ktrace.time(peak_idx:end), Iksum)
hold off
