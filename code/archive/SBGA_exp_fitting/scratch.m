clc
close all
clear variables

% read data
trace_data = table2array(readtable('4.5s-avg-ko.csv'));
t = trace_data(:, 1);

% voltage steps greater than 0 mv
yksum = trace_data(:, 7:end);
[~, num_volts] = size(yksum);

hold_idx = zeros(num_volts,1);
for i = 1:num_volts
    current_trc = yksum(:,i);
    [~, peak_idx] = max(current_trc);

    early_current_trc = current_trc(1:peak_idx);
    stable_val = min(early_current_trc);
    hold_idx(i) = find(early_current_trc == stable_val, 1, 'last');
end

%% test bi-exponential fitting
test_idx = 1;
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

%% test ikto
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

x = [30.0, 20.0, 33.5, 7.0, 0.4067];

volts = -50:10:50;
y = ikto(x, -70, volts(1), time_space, Ek);
plot(sim_t, y)

hold on
for i = 2:length(volts)
    y = ikto(x, -70, volts(i), time_space, Ek);
    plot(sim_t, y)
end
hold off

%% test ikslow
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

x = [22.5, 7.7, 45.2, 5.7, 1200, 0.16];

volts = -50:10:50;
y = ikslow(x, -70, volts(11), time_space, Ek);
plot(sim_t, y)

hold on
for i = 2:length(volts)
    y = ikslow(x, -70, volts(i), time_space, Ek);
    plot(sim_t, y)
end
hold off
