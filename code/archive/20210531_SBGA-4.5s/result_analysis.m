clc
close all
clear variables

% voltage clamp protocol
holding_p = -70; %mV
holding_t = 0.125*1000; %ms
P1_t = 4.5*1000; % ms
Ek = -91.1;
V = -50:10:50;
numv = length(V);

% load model fitting results - Ito WT
load('Ito_WT.mat')
to_param_wt = best_chroms;

[~, best_amps_idx] = min(best_amps);
[~, best_tau_idx] = min(best_taus);
if best_amps_idx == best_tau_idx
    best_toWt_idx = best_amps_idx;
else
    [~, best_toWt_idx] = min(best_amps+best_taus);
end

% load model fitting results - Ito KO
load('Ito_KO.mat')
to_param_ko = best_chroms;

[~, best_amps_idx] = min(best_amps);
[~, best_tau_idx] = min(best_taus);
if best_amps_idx == best_tau_idx
    best_toKo_idx = best_amps_idx;
else
    [~, best_toKo_idx] = min(best_amps+best_taus);
end

% load model fitting results - IKslow WT
load('IKslow_WT.mat')
kslow_param_wt = best_chroms;

[~, best_amps_idx] = min(best_amps);
[~, best_tau_idx] = min(best_taus);
if best_amps_idx == best_tau_idx
    best_KslowWt_idx = best_amps_idx;
else
    [~, best_KslowWt_idx] = min(best_amps+best_taus);
end

% load model fitting results - IKslow KO
load('IKslow_KO.mat')
kslow_param_ko = best_chroms;

[~, best_amps_idx] = min(best_amps);
[~, best_tau_idx] = min(best_taus);
if best_amps_idx == best_tau_idx
    best_KslowKo_idx = best_amps_idx;
else
    [~, best_KslowKo_idx] = min(best_amps+best_taus);
end

% clear temporary variables
clear best_amps best_amps_idx best_chroms best_tau_idx best_taus 

num_iters = 30;

%% compare with exponential function
% true amplitudes and taus
amps_wt = [24.93, 20.49];
taus_wt = [110.36, 1441.59];
amps_ko = [18.39, 4.70];
taus_ko = [111.86, 1874.92];

% time space for exponential function
t_holding = 0:0.1:holding_t;
t_P1 = 0:0.1:(P1_t-holding_t);
t_P1_shift = t_P1 + holding_t + 0.1;
t_exp = [t_holding, t_P1_shift];

% current for exponential function
holding_exp = zeros(1, length(t_holding)); 

Ito_wt_exp = exp_fn(t_P1, amps_wt(1), taus_wt(1));
Ito_wt_exp = [holding_exp, Ito_wt_exp];

Ito_ko_exp = exp_fn(t_P1, amps_ko(1), taus_ko(1));
Ito_ko_exp = [holding_exp, Ito_ko_exp];

IKslow_wt_exp = exp_fn(t_P1, amps_wt(2), taus_wt(2));
IKslow_wt_exp = [holding_exp, IKslow_wt_exp];

IKslow_ko_exp = exp_fn(t_P1, amps_ko(2), taus_ko(2));
IKslow_ko_exp = [holding_exp, IKslow_ko_exp];

figure(1)
subplot(1, 2, 1)
plot(t_exp, Ito_wt_exp, '--', 'LineWidth',2, 'Color','blue')
hold on
    for i = 1:num_iters
        [t_sim, ~, A1] = Ito(to_param_wt(i, :), holding_p, holding_t, 50, P1_t, Ek);
        plot(t_sim, A1(:, 5), 'LineWidth',1, 'Color',[0.1, 0.1, 0.1])
    end
hold off
axis tight
legend('Exponential fitting', 'Simulation Models')
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')

subplot(1, 2, 2)
plot(t_exp, IKslow_ko_exp, '--', 'LineWidth',2, 'Color','blue')
hold on
    for i = 1:num_iters
        [t_sim, ~, A1] = IKslow(kslow_param_ko(i, :), holding_p, holding_t, 50, P1_t, Ek);
        plot(t_sim, A1(:, 5), 'LineWidth',1, 'Color',[0.1, 0.1, 0.1, 0.1])
    end
hold off
axis tight
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')

figure(2)
subplot(1, 2, 1)
plot(t_exp, Ito_ko_exp, '--', 'LineWidth',2, 'Color','red')
hold on
    for i = 1:num_iters
        [t_sim, ~, A1] = Ito(to_param_ko(i, :), holding_p, holding_t, 50, P1_t, Ek);
        plot(t_sim, A1(:, 5), 'LineWidth',1, 'Color',[0.1, 0.1, 0.1, 0.1])
    end
hold off
axis tight
legend('Exponential Function', 'Simulation Models')
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')

subplot(1, 2, 2)
plot(t_exp, IKslow_ko_exp, '--', 'LineWidth',2, 'Color','red')
hold on
    for i = 1:num_iters
        [t_sim, ~, A1] = IKslow(kslow_param_ko(i, :), holding_p, holding_t, 50, P1_t, Ek);
        plot(t_sim, A1(:, 5), 'LineWidth',1, 'Color',[0.1, 0.1, 0.1, 0.1])
    end
hold off
axis tight
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')

%% prediction; generate traces with different voltage with the best
close all

% text string for legend
vtext = cell(1,numv);
for i=1:numv
    vtext{i} = strcat(num2str(V(numv-i+1)),' mV');
end
vtext = string(vtext);

% WT
figure('Color','w', 'Position',[100,100,480,200])
subplot(1, 2, 1)
hold on
for i=1:numv
    [t_to, ~, A_to] = Ito(to_param_wt(best_toWt_idx, :), holding_p, holding_t, V(numv-i+1), P1_t, Ek);
    plot(t_to, A_to(:, 5));
end
hold off
axis tight
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')
legend(vtext)
set(gca,'FontName','Arial','FontSize',11,'FontWeight','bold')

subplot(1, 2, 2)
hold on
for i = 1:numv
    [t_Kslow, ~, A_Kslow] = IKslow(kslow_param_wt(best_KslowWt_idx, :), holding_p, holding_t, V(numv-i+1), P1_t, Ek );
    plot(t_Kslow, A_Kslow(:, 5));
end
hold off
axis tight
ylabel('I_{Kslow} (pA/pF)')
xlabel('Time (ms)')
set(gca,'FontName','Arial','FontSize',11,'FontWeight','bold')

% Mgat1KO
figure('Color','w', 'Position',[600,100,480,200])
subplot(1, 2, 1)
hold on
for i = 1:numv
    [t_to, ~, A_to] = Ito(to_param_ko(best_toKo_idx, :), holding_p, holding_t, V(numv-i+1), P1_t, Ek);
    plot(t_to, A_to(:, 5));
end
hold off
axis tight
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')
legend(vtext)
set(gca,'FontName','Arial','FontSize',11,'FontWeight','bold')

subplot(1, 2, 2)
hold on
for i = 1:numv
    [t_Kslow, ~, A_Kslow] = IKslow(kslow_param_ko(best_KslowKo_idx, :), holding_p, holding_t, V(numv-i+1), P1_t, Ek );
    plot(t_Kslow, A_Kslow(:, 5));
end
hold off
axis tight
ylabel('I_{Kslow} (pA/pF)')
xlabel('Time (ms)')
set(gca,'FontName','Arial','FontSize',11,'FontWeight','bold')

%% prediction; SSA, SSI, densities, and time constants
% WT
to_peak_wt = zeros(num_iters, numv);
to_tau_wt = zeros(num_iters, numv);
to_ssa_wt = zeros(num_iters, numv);
to_ssi_wt = zeros(num_iters, numv);

kslow_peak_wt = zeros(num_iters, numv);
kslow_tau_wt = zeros(num_iters, numv);
kslow_ssa_wt = zeros(num_iters, numv);
kslow_ssi_wt = zeros(num_iters, numv);

% KO
to_peak_ko = zeros(num_iters, numv);
to_tau_ko = zeros(num_iters, numv);
to_ssa_ko = zeros(num_iters, numv);
to_ssi_ko = zeros(num_iters, numv);

kslow_peak_ko = zeros(num_iters, numv);
kslow_tau_ko = zeros(num_iters, numv);
kslow_ssa_ko = zeros(num_iters, numv);
kslow_ssi_ko = zeros(num_iters, numv);

for i = 1:num_iters
    for j = 1:numv
        % running simulation
        [to_twt, ~, to_Awt] = Ito(to_param_wt(i, :), holding_p, holding_t, V(j), P1_t, Ek);
        [kslow_twt, ~, kslow_Awt] = IKslow(kslow_param_wt(i, :), holding_p, holding_t, V(j), P1_t, Ek);

        [to_tko, ~, to_Ako] = Ito(to_param_ko(i, :), holding_p, holding_t, V(j), P1_t, Ek);
        [kslow_tko, ~, kslow_Ako] = IKslow(kslow_param_ko(i, :), holding_p, holding_t, V(j), P1_t, Ek);

        % peak (current densities)
        to_peak_wt(i, j) = max(to_Awt(:,5));
        kslow_peak_wt(i, j) = max(kslow_Awt(:,5));

        to_peak_ko(i, j) = max(to_Ako(:,5));
        kslow_peak_ko(i, j) = max(kslow_Ako(:,5));

        % inactivation time constant (tau)
        to_tau_wt(i, j) = 1/(to_Awt(end, 3)+to_Awt(end, 4));
        kslow_tau_wt(i, j) = kslow_Awt(end, 4);

        to_tau_ko(i, j) = 1/(to_Ako(end, 3)+to_Ako(end, 4));
        kslow_tau_ko(i, j) = kslow_Ako(end, 4);

        % SSA
        to_ssa_wt(i, j) = to_Awt(end,1)/(to_Awt(end,1)+to_Awt(end,2));
        kslow_ssa_wt(i, j) = kslow_Awt(end,1);
        
        to_ssa_ko(i, j) = to_Ako(end,1)/(to_Ako(end,1)+to_Ako(end,2));
        kslow_ssa_ko(i, j) = kslow_Ako(end,1);

        % SSI
        to_ssi_wt(i, j) = to_Awt(end, 3)/(to_Awt(end, 3)+to_Awt(end, 4));
        kslow_ssi_wt(i, j) = kslow_Awt(end, 2);

        to_ssi_ko(i, j) = to_Ako(end, 3)/(to_Ako(end, 3)+to_Ako(end, 4));
        kslow_ssi_ko(i, j) = kslow_Ako(end, 2);
    end
end

%% calculate standard error
% WT
to_peak_sem_wt = std(to_peak_wt)/sqrt(length(to_peak_wt));
to_tau_sem_wt = std(to_tau_wt)/sqrt(length(to_tau_wt));
to_ssa_sem_wt = std(to_ssa_wt)/sqrt(length(to_ssa_wt));
to_ssi_sem_wt = std(to_ssi_wt)/sqrt(length(to_ssi_wt));

kslow_peak_sem_wt = std(kslow_peak_wt)/sqrt(length(kslow_peak_wt));
kslow_tau_sem_wt = std(kslow_tau_wt)/sqrt(length(kslow_tau_wt));
kslow_ssa_sem_wt = std(kslow_ssa_wt)/sqrt(length(kslow_ssa_wt));
kslow_ssi_sem_wt = std(kslow_ssi_wt)/sqrt(length(kslow_ssi_wt));

% KO
to_peak_sem_ko = std(to_peak_ko)/sqrt(length(to_peak_ko));
to_tau_sem_ko = std(to_tau_ko)/sqrt(length(to_tau_ko));
to_ssa_sem_ko = std(to_ssa_ko)/sqrt(length(to_ssa_ko));
to_ssi_sem_ko = std(to_ssi_ko)/sqrt(length(to_ssi_ko));

kslow_peak_sem_ko = std(kslow_peak_ko)/sqrt(length(kslow_peak_ko));
kslow_tau_sem_ko = std(kslow_tau_ko)/sqrt(length(kslow_tau_ko));
kslow_ssa_sem_ko = std(kslow_ssa_ko)/sqrt(length(kslow_ssa_ko));
kslow_ssi_sem_ko = std(kslow_ssi_ko)/sqrt(length(kslow_ssi_ko));

%% current densities
figure(5)
subplot(1,2,1) % Ito
errorbar(V, mean(to_peak_wt, 1), to_peak_sem_wt, '--s', ...
'LineWidth',2, 'Color','blue')
hold on
    errorbar(V, mean(to_peak_ko, 1), to_peak_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
xlabel('Voltage (mV)')
ylabel('A_{to} (pA/pF)')
legend('WT', 'Mgat1KO')

subplot(1,2,2) % IKslow
errorbar(V, mean(kslow_peak_wt, 1), kslow_peak_sem_wt, '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V, mean(kslow_peak_ko, 1), kslow_peak_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
xlabel('Voltage (mV)')
ylabel('A_{Kslow} (pA/pF)')

%% SSA
figure(6)
subplot(1,2,1) % Ito
errorbar(V, to_ssa_wt(1,:), to_ssa_sem_wt, '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V, to_ssa_ko(1,:), to_ssa_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
ylim([0, 1])
xlabel('Voltage (mV)')
ylabel('SSA_{to}')
legend('WT', 'Mgat1KO')

subplot(1,2,2) % IKslow
errorbar(V, mean(kslow_ssa_wt, 1), kslow_ssa_sem_wt, '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V, mean(kslow_ssa_ko, 1), kslow_ssa_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
ylim([0, 1])
xlabel('Voltage (mV)')
ylabel('SSA_{Kslow}')
legend('WT', 'Mgat1KO')

%% inactivationtime constants
figure(7)
mean_to_tau_wt = mean(to_tau_wt, 1);
mean_to_tau_ko = mean(to_tau_ko, 1);
subplot(1, 2, 1) % Ito
errorbar(V(3:end), mean_to_tau_wt(3:end), to_tau_sem_wt(3:end), '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V(3:end), mean_to_tau_ko(3:end), to_tau_sem_ko(3:end), '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
xlabel('Voltage (mV)')
ylabel('\tau_{to} (ms)')
legend('WT', 'Mgat1KO')

subplot(1, 2, 2) % IKslow
mean_kslow_tau_wt = mean(kslow_tau_wt, 1);
mean_kslow_tau_ko = mean(kslow_tau_ko, 1);
errorbar(V, mean(kslow_tau_wt, 1), kslow_tau_sem_wt, '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V, mean(kslow_tau_ko, 1), kslow_tau_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
xlabel('Voltage (mV)')
ylabel('\tau_{to} (ms)')

%% SSI
figure(8)
subplot(1, 2, 1) % Ito
errorbar(V, mean(to_ssi_wt, 1), to_ssi_sem_wt, '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V, mean(to_ssi_ko, 1), to_ssi_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
ylim([0, 1])
xlabel('Voltage (mV)')
ylabel('SSI_{to}')
legend('WT', 'Mgat1KO')

subplot(1, 2, 2) % IKslow
errorbar(V, mean(kslow_ssi_wt, 1), kslow_ssi_sem_wt, '--s', ...
'Color','blue', 'LineWidth',2)
hold on
    errorbar(V, mean(kslow_ssi_ko, 1), kslow_ssi_sem_ko, '-o', ...
    'Color','red', 'LineWidth',2)
hold off
axis tight
grid on
xticks(V)
ylim([0, 1])
xlabel('Voltage (mV)')
ylabel('SSI_{Kslow}')

%% bar plot for 4.5-sec protocol data
clc
close all
clear variables

wt = readtable('Ik-4.5s-KO.csv');
ko = readtable('Ik-4.5s-WT.csv');

% mean WT
to_wt = mean(wt.peak_Ito);
kslow_wt = mean(wt.peak_IKslow);
ss_wt = mean(wt.Iss);
tot_wt = mean(wt.tau_Ito);
kslowt_wt = mean(wt.tau_IKslow);

% error WT
to_eb_wt = std(wt.peak_Ito)/sqrt(length(wt.peak_Ito)); 
kslow_eb_wt = std(wt.peak_IKslow)/sqrt(length(wt.peak_IKslow));
ss_eb_wt = std(wt.Iss)/sqrt(length(wt.Iss));
tot_eb_wt = std(wt.tau_Ito)/sqrt(length(wt.tau_Ito)); 
kslowt_eb_wt = std(wt.tau_IKslow)/sqrt(length(wt.tau_IKslow)); 

% mean KO
to_ko = mean(ko.peak_Ito);
kslow_ko = mean(ko.peak_IKslow);
ss_ko = mean(ko.Iss);
tot_ko = mean(ko.tau_Ito);
kslowt_ko = mean(ko.tau_IKslow);

% error KO
to_eb_ko = std(ko.peak_Ito)/sqrt(length(ko.peak_Ito)); 
kslow_eb_ko = std(ko.peak_IKslow)/sqrt(length(ko.peak_IKslow)); 
ss_eb_ko = std(ko.Iss)/sqrt(length(ko.Iss));
tot_eb_ko = std(ko.tau_Ito)/sqrt(length(ko.tau_Ito)); 
kslowt_eb_ko = std(ko.tau_IKslow)/sqrt(length(ko.tau_IKslow)); 

vals1 = zeros(3, 2);
errs1 = zeros(3, 2);

vals2 = zeros(2, 2);
errs2 = zeros(2, 2);

vals1(1, :) = [to_wt, to_ko];
vals1(2, :) = [kslow_wt, kslow_ko];
vals1(3, :) = [ss_wt, ss_ko];
errs1(1, :) = [to_eb_wt, to_eb_ko];
errs1(2, :) = [kslow_eb_wt, kslow_eb_ko];
errs1(3, :) = [ss_eb_wt, ss_eb_ko];

vals2(1, :) = [tot_wt, tot_ko];
vals2(2, :) = [kslowt_wt, kslowt_ko];
errs2(1, :) = [tot_eb_wt, tot_eb_ko];
errs2(2, :) = [kslowt_eb_wt, kslowt_eb_ko];

% peaks
figure(1)
b = bar(vals1, 'grouped');
num_bars = 2;
x = [];
hold on
    for i = 1:num_bars
        x = [x; b(i).XEndPoints];
    end
errorbar(x', vals1, errs1, 'k', 'LineStyle','none')
hold off
xticklabels({'I_{to}','I_{Kslow}','I_{SS}'})
ylabel('Peak (pA/pF)')
legend('WT','Mgat1KO')

% taus
figure(2)
b = bar(vals2, 'grouped');
num_bars = 2;
x = [];
hold on
    for i = 1:num_bars
        x = [x; b(i).XEndPoints];
    end
errorbar(x', vals2, errs2, 'k', 'LineStyle','none')
hold off
xticklabels({'\tau_{to}','\tau_{Kslow}'})
ylabel('Time Constant (ms)')

%% fitted parameters mean and SEM
disp(mean(to_param_wt))
disp(std(to_param_wt)/sqrt(length(to_param_wt)))

disp(mean(to_param_ko))
disp(std(to_param_ko)/sqrt(length(to_param_ko)))

disp(mean(kslow_param_wt))
disp(std(kslow_param_wt)/sqrt(length(kslow_param_wt)))

disp(mean(kslow_param_ko))
disp(std(kslow_param_ko)/sqrt(length(kslow_param_ko)))

%% prediction: peaks at different voltages
exp_ksum = table2array(readtable('./4.5s-avg-wt.csv'));
exp_t = exp_ksum(:,1);
exp_yksum = exp_ksum(:,7:end);

volts = 0:10:50;
exp_peaks = zeros(length(volts), 1);
sim_peaks = zeros(length(volts), num_iters);
for i = 1:length(volts)
    exp_peaks(i) = max(exp_yksum(:, i));
    for j = 1:num_iters
       [t_to, ~, A1] = Ito(to_param_wt(j, :), holding_p, holding_t, volts(i), P1_t, Ek);
       [t_kslow, ~, A2] = IKslow(kslow_param_ko(j, :), holding_p, holding_t, volts(i), P1_t, Ek);
       
       peak_to = max(A1(:, 5));
       peak_kslow = max(A2(:, 5));
       
       sim_peaks(i, j) = peak_to + peak_kslow;
    end
end
