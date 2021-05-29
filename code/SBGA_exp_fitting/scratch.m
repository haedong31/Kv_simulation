clc
close all
clear variables

% read data
trace_data = table2array(readtable('4.5s-avg-wt.csv'));
t = trace_data(:, 1);

% voltage steps greater than 0 mv
yksum = trace_data(:, 7:end);
[~, num_volts] = size(yksum);

%% test bi-exponential fitting
test_idx = 6;
testy = yksum(:, test_idx);

[~, peak_idx] = max(testy);
t_trunc = t(peak_idx:end);
t_trunc = t_trunc - t_trunc(1);

% bi-exponential fitting
[amp, tau] = bi_exp_fit(t, testy);

bi_exp = @(amp, tau) amp(1).*exp(-t_trunc./tau(1)) + amp(2).*exp(-t_trunc./tau(2)) + amp(3);
yhat = bi_exp(amp, tau);

% visualization
plot(t_trunc, testy(peak_idx:end));
hold on
plot(t_trunc, yhat)
hold off
legend('experimental','fitting')

