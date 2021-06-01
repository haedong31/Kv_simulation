clc
close all
clear variables

num_iters = 30;

% load model fitting results and select best one
load('calib_result_wt.mat')

load('calib_result_ko.mat')

% ikto wt
[~, best_amps_idx] = min(to_amps_wt);
[~, best_tau_idx] = min(to_taus_wt);
if best_amps_idx == best_tau_idx
    best_toWt_idx = best_amps_idx;
else
    [~, best_toWt_idx] = min(to_amps_wt+to_taus_wt);
end

% ikslow wt
[~, best_amps_idx] = min(kslow_amps_wt);
[~, best_tau_idx] = min(kslow_taus_wt);
if best_amps_idx == best_tau_idx
    best_KslowWt_idx = best_amps_idx;
else
    [~, best_KslowWt_idx] = min(kslow_amps_wt+kslow_taus_wt);
end

% ikto ko
[~, best_amps_idx] = min(to_amps_ko);
[~, best_tau_idx] = min(to_taus_ko);
if best_amps_idx == best_tau_idx
    best_toKo_idx = best_amps_idx;
else
    [~, best_toKo_idx] = min(to_amps_ko+to_taus_ko);
end

% ikslow ko
[~, best_amps_idx] = min(kslow_amps_ko);
[~, best_tau_idx] = min(kslow_taus_ko);
if best_amps_idx == best_tau_idx
    best_KslowKo_idx = best_amps_idx;
else
    [~, best_KslowKo_idx] = min(kslow_amps_ko+kslow_taus_ko);
end

% clear temporary variables
clear best_amps_idx best_tau_idx 

%% compare with experimental data
volts = 0:10:50;
num_volts = length(volts);

trace_data_wt = table2array(readtable('4.5s-avg-wt.csv'));
trace_data_ko = table2array(readtable('4.5s-avg-ko.csv'));

t_wt = trace_data_wt(:, 1);
t_ko = trace_data_ko(:, 1);

yksum_wt = trace_data_wt(:, 7:end);
yksum_ko = trace_data_ko(:, 7:end);

% estimate holding time
[~, ideal_hold_idx_wt] = min(abs(t_wt-120));
[~, ideal_hold_idx_ko] = min(abs(t_ko-120));

init_stable_val_wt = sqrt(var(yksum_wt(ideal_hold_idx_wt,:)));
init_stable_val_ko = sqrt(var(yksum_wt(ideal_hold_idx_ko,:)));

hold_idx_wt = zeros(num_volts, 1);
hold_idx_ko = zeros(num_volts, 1);

for i = 1:num_volts
    current_trace_wt = yksum_wt(:, i);

    counter = 1;
    statbility_est = abs(current_trace_wt(ideal_hold_idx_wt+counter) - current_trace_wt(ideal_hold_idx_wt-counter));
    
    if statbility_est > 10*init_stable_val_wt
        hold_idx_wt(i) = ideal_hold_idx_wt;
        continue
    end

    while true
        % update counter and stable value
        counter = counter + 1;
        stable_val = statbility_est;
        
        % check stability
        statbility_est = abs(current_trace_wt(ideal_hold_idx_wt+counter) - current_trace_wt(ideal_hold_idx_wt-counter));
        if statbility_est > 10*stable_val
            hold_idx_wt(i) = ideal_hold_idx_wt + (counter-1);
            break
        end
    end
end

for i = 1:num_volts
    current_trace_ko = yksum_ko(:, i);

    counter = 1;
    statbility_est = abs(current_trace_ko(ideal_hold_idx_ko+counter) - current_trace_ko(ideal_hold_idx_ko-counter));
    
    if statbility_est > 10*init_stable_val_ko
        hold_idx_ko(i) = ideal_hold_idx_ko;
        continue
    end

    while true
        % update counter and stable value
        counter = counter + 1;
        stable_val = statbility_est;
        
        % check stability
        statbility_est = abs(current_trace_ko(ideal_hold_idx_ko+counter) - current_trace_ko(ideal_hold_idx_ko-counter));
        if statbility_est > 10*stable_val
            hold_idx_ko(i) = ideal_hold_idx_ko + (counter-1);
            break
        end
    end
end

%% visual inspection
Ek = -91.1;

par_kto = to_param_wt(best_toWt_idx, :);
par_kslow = kslow_param_wt(best_KslowWt_idx, :);

t = t_wt;

figure(1)
volt_idx = 1;
hold_pt = t(hold_idx_wt(volt_idx));
end_pt = t(end);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx_wt(volt_idx));
time_space{3} = t(hold_idx_wt(volt_idx)+1:end) - t(hold_idx_wt(volt_idx)+1);

y1 = ikto(par_kto, -70, volts(1), time_space, Ek);
y2 = ikslow(par_kslow, -70, volts(1), time_space, Ek);
y3 = zeros(length(t), 1);
y3(hold_idx_wt(volt_idx)+1:end) = amp_wt(3,1);

plot(t, yksum_wt(:, 1))
hold on
    plot(t, y1+y2+y3)
hold off
legend('Experiemntal', 'Simulated')

figure(2)
y1 = ikto(par_kto, -70, volts(2), time_space, Ek);
y2 = ikslow(par_kslow, -70, volts(2), time_space, Ek);
y3 = zeros(length(sim_t), 1);
y3((hold_len + 1):end_len) = amp(3,2);

plot(t, yksum(:, 2))
hold on
    plot(sim_t, y1+y2+y3)
hold off
legend('Experiemntal', 'Simulated')

figure(3)
y1 = ikto(par_kto, -70, volts(3), time_space, Ek);
y2 = ikslow(par_kslow, -70, volts(3), time_space, Ek);
y3 = zeros(length(sim_t), 1);
y3((hold_len + 1):end_len) = amp(3,3);

plot(t, yksum(:, 3))
hold on
    plot(sim_t, y1+y2+y3)
hold off
legend('Experiemntal', 'Simulated')

figure(4)
y1 = ikto(par_kto, -70, volts(4), time_space, Ek);
y2 = ikslow(par_kslow, -70, volts(4), time_space, Ek);
y3 = zeros(length(sim_t), 1);
y3((hold_len + 1):end_len) = amp(3,4);

plot(t, yksum(:, 4))
hold on
    plot(sim_t, y1+y2+y3)
hold off
legend('Experiemntal', 'Simulated')

figure(5)
y1 = ikto(par_kto, -70, volts(5), time_space, Ek);
y2 = ikslow(par_kslow, -70, volts(5), time_space, Ek);
y3 = zeros(length(sim_t), 1);
y3((hold_len + 1):end_len) = amp(3,5);

plot(t, yksum(:, 5))
hold on
    plot(sim_t, y1+y2+y3)
hold off
legend('Experiemntal', 'Simulated')

figure(6)
volt_idx = 6;
hold_pt = t(hold_idx_wt(volt_idx));
end_pt = t(end);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx_wt(volt_idx));
time_space{3} = t(hold_idx_wt(volt_idx)+1:end) - t(hold_idx_wt(volt_idx)+1);

y1 = ikto(par_kto, -70, volts(volt_idx), time_space, Ek);
y2 = ikslow(par_kslow, -70, volts(volt_idx), time_space, Ek);
y3 = zeros(length(t), 1);
y3(hold_idx_wt(volt_idx)+1:end) = amp_wt(3,volt_idx);

plot(t, yksum_wt(:, volt_idx))
hold on
    plot(t, y1+y2+y3)
hold off
legend('Experiemntal', 'Simulated')

%% compare with exponential function
% protocol
end_len = 4.5*1000;
hold_len = 0.125*1000;

hold_t = 0:hold_len;
pulse_t = (hold_len + 1):end_len;
pulse_t_adj = pulse_t - pulse_t(1);
sim_t = [hold_t, pulse_t];

time_space = cell(1,3);
time_space{1} = sim_t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

Ek = -91.1;

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
figure(3) % WT
subplot(1, 2, 1)
hold on
    for i = 1:numv
        [t_to, ~, A_to] = Ito(to_param_wt(best_toWt_idx, :), holding_p, holding_t, V(i), P1_t, Ek);
        plot(t_to, A_to(:, 5));
    end
hold off
axis tight
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')

subplot(1, 2, 2)
hold on
    for i = 1:numv
        [t_Kslow, ~, A_Kslow] = IKslow(kslow_param_wt(best_KslowWt_idx, :), holding_p, holding_t, V(i), P1_t, Ek );
        plot(t_Kslow, A_Kslow(:, 5));
    end
hold off
axis tight
ylabel('I_{Kslow} (pA/pF)')
xlabel('Time (ms)')

figure(4) % KO
subplot(1, 2, 1)
hold on
    for i = 1:numv
        [t_to, ~, A_to] = Ito(to_param_ko(best_toKo_idx, :), holding_p, holding_t, V(i), P1_t, Ek);
        plot(t_to, A_to(:, 5));
    end
hold off
axis tight
ylabel('I_{to} (pA/pF)')
xlabel('Time (ms)')

subplot(1, 2, 2)
hold on
    for i = 1:numv
        [t_Kslow, ~, A_Kslow] = IKslow(kslow_param_ko(best_KslowKo_idx, :), holding_p, holding_t, V(i), P1_t, Ek );
        plot(t_Kslow, A_Kslow(:, 5));
    end
hold off
axis tight
ylabel('I_{Kslow} (pA/pF)')
xlabel('Time (ms)')

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
