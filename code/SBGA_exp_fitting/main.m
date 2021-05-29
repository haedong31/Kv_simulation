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

tol = [0.1, 1];
N0 = 30;
N1 = 6;
N2 = 4;
par = ikto_calibration(amp_kto, tau_kto, tol, N0, N1, N2);
