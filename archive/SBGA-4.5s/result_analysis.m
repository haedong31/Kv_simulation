%% verification
clc
close all
clear variables

% load model fitting results 
num_iters = 30;
load('calib_result_wt.mat')
load('calib_result_ko.mat')

% read data
trace_data_wt = table2array(readtable('4.5s-avg-wt.csv'));
trace_data_ko = table2array(readtable('4.5s-avg-ko.csv'));

t_wt = trace_data_wt(:, 1);
t_ko = trace_data_ko(:, 1);

% protocol and estimate holding time index
hold_volt = -70;
Ek = -91.1;
ideal_hold_time = 0.125*1000;
[~, hold_idx_wt] = min(abs(t_wt - ideal_hold_time));
[~, hold_idx_ko] = min(abs(t_ko - ideal_hold_time));

% exclude data with voltage steps less than 0 mv
volts = 0:10:50;
num_volts = length(volts);
yksum_wt = trace_data_wt(:, 7:end);
yksum_ko = trace_data_ko(:, 7:end);

% simulation time space
pulse_t = t_wt((hold_idx_wt+1):end);
pulse_t_adj = pulse_t - pulse_t(1);

time_space_wt = cell(1,3);
time_space_wt{1} = t_wt;
time_space_wt{2} = t_wt(1:hold_idx_wt);
time_space_wt{3} = pulse_t_adj;

pulse_t = t_ko((hold_idx_ko+1):end);
pulse_t_adj = pulse_t - pulse_t(1);

time_space_ko = cell(1,3);
time_space_ko{1} = t_ko;
time_space_ko{2} = t_ko(1:hold_idx_ko);
time_space_ko{3} = pulse_t_adj;

tau_diff_wt = zeros(num_iters, num_volts);
for i = 1:num_iters
    for j = 1:num_volts
        volt_idx = j;
        volt = volts(volt_idx);
        
        par_kto = to_param_wt(i, :);
        par_kslow = kslow_param_wt(i, :);
        
        ykto = ikto(par_kto, hold_volt, volt, time_space_wt, Ek);
        ykslow = ikslow(par_kslow, hold_volt, volt, time_space_wt, Ek);
        
        yksum = yksum_wt(:, volt_idx);
        yksum_hat = ykto + ykslow;

        [peak, peak_idx] = max(yksum);
        [peak_hat, peak_hat_idx] = max(yksum_hat);

        yksum_trunc = yksum(peak_idx:end);
        t_wt_trunc = t_wt(peak_idx:end);
        t_wt_trunc = t_wt_trunc - t_wt_trunc(1);
        [~, tau_idx] = min(abs(peak*exp(-1) - yksum_trunc));
        tau = t_wt_trunc(tau_idx);

        yksum_hat_trunc = yksum_hat(peak_hat_idx:end);
        t_wt_trunc = t_wt(peak_hat_idx:end);
        t_wt_trunc = t_wt_trunc - t_wt_trunc(1);
        [~, tau_hat_idx] = min(abs(peak_hat*exp(-1) - yksum_hat_trunc));
        tau_hat = t_wt_trunc(tau_hat_idx);

        tau_diff_wt(i, j) = abs(tau - tau_hat);
    end
end

weight = exp(0:(1/(num_volts-1)):1);
weight = weight/sum(weight);
[~, best_model_idx_wt] = min(mean(weight.*tau_diff_wt, 2));

volt_idx = 2;
volt = volts(volt_idx);

par_kto = to_param_wt(best_model_idx_wt, :);
par_kslow = kslow_param_wt(best_model_idx_wt, :);

yksum = yksum_wt(:, volt_idx);
ykto = ikto(par_kto, hold_volt, volt, time_space_wt, Ek);
ykslow = ikslow(par_kslow, hold_volt, volt, time_space_wt, Ek);
yksum_hat = ykto + ykslow ;

figure(1)
plot(t_wt, yksum)
hold on
plot(t_wt, yksum_hat)
hold off
legend('Experimental','Simulatied')

tau_diff_ko = zeros(num_iters, num_volts);
for i = 1:num_iters
    for j = 1:num_volts
        volt_idx = j;
        volt = volts(volt_idx);
        
        par_kto = to_param_ko(i, :);
        par_kslow = kslow_param_ko(i, :);
        
        ykto = ikto(par_kto, hold_volt, volt, time_space_ko, Ek);
        ykslow = ikslow(par_kslow, hold_volt, volt, time_space_ko, Ek);
        
        yksum = yksum_ko(:, volt_idx);
        yksum_hat = ykto + ykslow;

        [peak, peak_idx] = max(yksum);
        [peak_hat, peak_hat_idx] = max(yksum_hat);

        yksum_trunc = yksum(peak_idx:end);
        t_ko_trunc = t_ko(peak_idx:end);
        t_ko_trunc = t_ko_trunc - t_ko_trunc(1);
        [~, tau_idx] = min(abs(peak*exp(-1) - yksum_trunc));
        tau = t_ko_trunc(tau_idx);

        yksum_hat_trunc = yksum_hat(peak_hat_idx:end);
        t_ko_trunc = t_ko(peak_hat_idx:end);
        t_ko_trunc = t_ko_trunc - t_ko_trunc(1);
        [~, tau_hat_idx] = min(abs(peak_hat*exp(-1) - yksum_hat_trunc));
        tau_hat = t_ko_trunc(tau_hat_idx);

        tau_diff_ko(i, j) = abs(tau - tau_hat);
    end
end

weight = exp(0:(1/(num_volts-1)):1);
weight = weight/sum(weight);
[~, best_model_idx_ko] = min(mean(weight.*tau_diff_ko, 2));

volt_idx = 4;
volt = volts(volt_idx);

par_kto = to_param_ko(best_model_idx_ko, :);
par_kslow = kslow_param_ko(best_model_idx_ko, :);

yksum = yksum_ko(:, volt_idx);
ykto = ikto(par_kto, hold_volt, volt, time_space_ko, Ek);
ykslow = ikslow(par_kslow, hold_volt, volt, time_space_ko, Ek);
yksum_hat = ykto + ykslow ;

figure(2)
plot(t_ko, yksum)
hold on
plot(t_ko, yksum_hat)
hold off
legend('Experimental','Simulatied')
