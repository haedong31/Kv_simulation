clc
close all
clear variables


%% data
% protocol
holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
P2 = -70; % mV
P2_t = P1_t; % ms

% raw trace & normalization with capacitance
load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
% load('ds_Ktrace_wt.mat')
% ds_Ktrace = ds_Ktrace_wt;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};

K_data = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx');
% K_data = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx');

Iss = K_data.IssFF;
% Iss = K_data.Iss;
Iss = nanmean(Iss);

Ito = K_data.A3FF;
% Ito = K_data.A3;
Ito = nanmean(Ito);

tau_to = K_data.Tau3FF;
% tau_to = K_data.Tau3;
tau_to = nanmean(tau_to);

IKslow1 = K_data.A2FF;
% IKslow1 = K_data.A2;
IKslow1 = nanmean(IKslow1);

tau1 = K_data.Tau2FF;
% tau1 = K_data.Tau2;
tau1 = nanmean(tau1);

IKslow2 = K_data.A1FF;
% IKslow2 = K_data.A1;
IKslow2 = nanmean(IKslow2);

tau2 = K_data.Tau1FF;
% tau2 = K_data.Tau1;
tau2 = nanmean(tau2);

cap = K_data.CapFF;
% cap = K_data.Cap;
cap = nanmean(cap);

ds_Ktrace.I = ds_Ktrace.I ./ cap;


%% run simulation
[t, ~, A, ~] = Kv(holding_p, holding_t, P1, P1_t, P2, P2_t);

% IKsum without Iss
Ito_trc = A(:,13);
IKslow1_trc = A(:,14);
IKslow2_trc = A(:,15);
IKsum = Ito_trc + IKslow1_trc + IKslow2_trc;

% add the constant term (Iss) after the peak
[peak, peak_idx] = max(IKsum);
IKsum(peak_idx:end) = IKsum(peak_idx:end) + Iss;

figure(1)
plot(t, IKsum)
hold on
plot(ds_Ktrace.time, ds_Ktrace.I, 'LineWidth',2)
hold off
legend('IKsum', 'Raw Trace')

% align the two traces and downsampling IKsum
IKsum = alignsignals(IKsum, ds_Ktrace.I);
IKsum_tbl = table(t, IKsum);
sample_rate = floor(length(t)/length(ds_Ktrace.time));

if sample_rate > 0
    ds_IKsum = downsample(IKsum_tbl, sample_rate);
else
    ds_IKsum = Iksum_tbl;
end

figure(2)
plot(ds_IKsum.t, ds_IKsum.IKsum)
hold on
plot(ds_Ktrace.time, ds_Ktrace.I, 'LineWidth',2)
hold off
title('Aligned')
legend('IKsum', 'Raw Trace')


%% calculate objective function
trace_sim = dtw(ds_IKsum.IKsum, ds_Ktrace.I);

Ito_hat = max(Ito_trc);
IKslow1_hat = max(IKslow1_trc);
IKslow2_hat = max(IKslow2_trc);

amp_to = abs(Ito - max(Ito_trc));
amp_slow1 = abs(IKslow1 - max(IKslow1_trc));
amp_slow2 = abs(IKslow2 - max(IKslow2_trc));

[~, tau_to_idx] = min(abs(Ito_hat*exp(-1) - Ito_trc));
[~, tau1_idx] = min(abs(IKslow1_hat*exp(-1) - IKslow1_trc));
[~, tau2_idx] = min(abs(IKslow2_hat*exp(-1) - IKslow2_trc));

tau_to = abs(tau_to - t(tau_to_idx));
tau_slow1 = abs(tau1 - t(tau1_idx));
tau_slow2 = abs(tau2 - t(tau2_idx));

z = trace_sim + amp_to + amp_slow1 + amp_slow2 + tau_to + tau_slow1 + tau_slow2;
