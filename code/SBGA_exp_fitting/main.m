clc
close all
clear variables

% read data
trace_data = table2array(readtable('4.5s-avg-wt.csv'));
t = trace_data(:, 1);

% exclude data with voltage steps less than 0 mv
yksum = trace_data(:, 7:end);
volts = 0:10:50;

% apply bi-exponential fitting
num_volts = length(volts);
amp = zeros(3, num_volts);
tau = zeros(2, num_volts);

for i = 1:num_volts
    [amp_running, tau_running] = bi_exp_fit(t, yksum(:, i));
    amp(:, i) = amp_running;
    tau(:, i) = tau_running;
end

% classify amp and tau of IKto and IKslow
amp_kto = zeros(num_volts, 1);
amp_kslow = zeros(num_volts, 1);
tau_kto = zeros(num_volts, 1);
tau_kslow = zeros(num_volts, 1);

for i = 1:num_volts
    tau1 = tau(1, i);
    tau2 = tau(2, i);

    if tau1 < tau2
        amp_kto(i) = amp(1, i);
        amp_kslow(i) = amp(2, i);

        tau_kto(i) = tau(1, i);
        tau_kslow(i) = tau(2, i);
    else
        amp_kto(i) = amp(2, i);
        amp_kslow(i) = amp(1, i);

        tau_kto(i) = tau(2, i);
        tau_kslow(i) = tau(1, i);
    end
end

%% calibration
% IKto
num_iters = 30;
N1 = 10;
N2 = 9;

% run calibration
to_param_wt = zeros(num_iters, 5);
to_amps_wt = zeros(num_iters, 1);
to_taus_wt = zeros(num_iters, 1);

kslow_param_wt = zeros(num_iters, 6);
kslow_amps_wt = zeros(num_iters, 1);
kslow_taus_wt = zeros(num_iters, 1);

for i = 1:num_iters
    fprintf('[%i/%i] \n', i, num_iters)

    [par, amp_diff, tau_diff] = ikto_calibration(amp_kto, tau_kto, volts, N1, N2, false);
    to_param_wt(i, :) = par;
    to_amps_wt(i) = amp_diff;
    to_taus_wt(i) = tau_diff;

    [par, amp_diff, tau_diff] = ikslow_calibration(amp_kslow, tau_kslow, volts, N1, N2, false);
    kslow_param_wt(i, :) = par;
    kslow_amps_wt(i) = amp_diff;
    kslow_taus_wt(i) = tau_diff;
end

save('calib_result_wt.mat', 'to_param_wt', 'to_amps_wt', 'to_taus_wt', 'kslow_param_wt', 'kslow_amps_wt', 'kslow_taus_wt')
