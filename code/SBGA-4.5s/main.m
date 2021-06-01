clc
close all
clear variables

% read data
trace_data = table2array(readtable('4.5s-avg-wt.csv'));
t = trace_data(:, 1);

% exclude data with voltage steps less than 0 mv
yksum = trace_data(:, 7:end);

% apply bi-exponential fitting
[amps, taus] = bi_exp_fit(t, yksum(:, end));

% classify amp and tau of IKto and IKslow
tau1 = taus(1);
tau2 = taus(2);

if tau1 < tau2
    amp_kto = amps(1);
    amp_kslow = amps(2);

    tau_kto = taus(1);
    tau_kslow = taus(2);
else
    amp_kto = amps(2);
    amp_kslow = amps(1);

    tau_kto = taus(2);
    tau_kslow = taus(1);
end

% calibration
% IKto
num_iters = 30;

input_volt = 50;
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

    [par, amp_diff, tau_diff] = ikto_calibration(amp_kto, tau_kto, input_volt, N1, N2);
    to_param_wt(i, :) = par;
    to_amps_wt(i) = amp_diff;
    to_taus_wt(i) = tau_diff;

    [par, amp_diff, tau_diff] = ikslow_calibration(amp_kslow, tau_kslow, input_volt, N1, N2);
    kslow_param_wt(i, :) = par;
    kslow_amps_wt(i) = amp_diff;
    kslow_taus_wt(i) = tau_diff;
end

amp_wt = amps;
save('calib_result_wt.mat', 'to_param_wt', 'to_amps_wt', 'to_taus_wt', 'kslow_param_wt', 'kslow_amps_wt', 'kslow_taus_wt', 'amp_wt')

%% same calibration routine for Mgat1KO
clc
close all
clear variables

% read data
trace_data = table2array(readtable('4.5s-avg-ko.csv'));
t = trace_data(:, 1);

% exclude data with voltage steps less than 0 mv
yksum = trace_data(:, 7:end);

% apply bi-exponential fitting
[amps, taus] = bi_exp_fit(t, yksum(:, end));

% classify amp and tau of IKto and IKslow
tau1 = taus(1);
tau2 = taus(2);

if tau1 < tau2
    amp_kto = amps(1);
    amp_kslow = amps(2);

    tau_kto = taus(1);
    tau_kslow = taus(2);
else
    amp_kto = amps(2);
    amp_kslow = amps(1);

    tau_kto = taus(2);
    tau_kslow = taus(1);
end

% calibration
% IKto
num_iters = 30;

input_volt = 50;
N1 = 10;
N2 = 9;

% run calibration
to_param_ko = zeros(num_iters, 5);
to_amps_ko = zeros(num_iters, 1);
to_taus_ko = zeros(num_iters, 1);

kslow_param_ko = zeros(num_iters, 6);
kslow_amps_ko = zeros(num_iters, 1);
kslow_taus_ko = zeros(num_iters, 1);

for i = 1:num_iters
    fprintf('[%i/%i] \n', i, num_iters)

    [par, amp_diff, tau_diff] = ikto_calibration(amp_kto, tau_kto, input_volt, N1, N2);
    to_param_ko(i, :) = par;
    to_amps_ko(i) = amp_diff;
    to_taus_ko(i) = tau_diff;

    [par, amp_diff, tau_diff] = ikslow_calibration(amp_kslow, tau_kslow, input_volt, N1, N2);
    kslow_param_ko(i, :) = par;
    kslow_amps_ko(i) = amp_diff;
    kslow_taus_ko(i) = tau_diff;
end

amp_ko = amps;
save('calib_result_ko.mat', 'to_param_ko', 'to_amps_ko', 'to_taus_ko', 'kslow_param_ko', 'kslow_amps_ko', 'kslow_taus_ko', 'amp_ko')
