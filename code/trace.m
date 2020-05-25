clc
close all
clear variables


%% manipulate the average raw trace
Ktrace = readtable('k_trace.csv');
Ktrace_ko = Ktrace(:, [1,2]);
Ktrace_wt = Ktrace(:, [1,3]);

plot(Ktrace_ko.time, Ktrace_ko.KO, 'LineWidth',2)
hold on
plot(Ktrace_wt.time, Ktrace_wt.WT, 'LineWidth',2)
hold off


%% downsample the raw trace
% cut off the raw traces
sub_Ktrace_ko = Ktrace_ko(1:500001,:);
sub_Ktrace_wt = Ktrace_wt(1:500001,:);

figure(1)
plot(sub_Ktrace_ko.time, sub_Ktrace_ko.KO, 'LineWidth',2)
hold on
plot(sub_Ktrace_wt.time, sub_Ktrace_wt.WT, 'LineWidth',2)
hold off

[peak_ko, peak_ko_idx] = max(sub_Ktrace_ko.KO);
[peak_wt, peak_wt_idx] = max(sub_Ktrace_wt.WT);

ds_Ktrace_ko = downsample(sub_Ktrace_ko, 98);
ds_Ktrace_wt = downsample(sub_Ktrace_wt, 98);

ds_Ktrace_ko.KO(ds_Ktrace_ko.KO < 0) = 0;
ds_Ktrace_wt.WT(ds_Ktrace_wt.WT < 0) = 0;

figure(2)
plot(ds_Ktrace_ko.time, ds_Ktrace_ko.KO, 'LineWidth',2)
hold on
plot(ds_Ktrace_wt.time, ds_Ktrace_wt.WT, 'LineWidth',2)
hold off
axis tight


%% amplutes and taus for IKslow1 and 2
Kko = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx');
Kwt = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx');

IKslow1_ko = Kko.A2FF;
IKslow1_ko = nanmean(IKslow1_ko);
IKslow1_wt = Kwt.A2;
IKslow1_wt = nanmean(IKslow1_wt);

IKslow2_ko = Kko.A1FF;
IKslow2_ko = nanmean(IKslow2_ko);
IKslow2_wt = Kwt.A1;
IKslow2_wt = nanmean(IKslow2_wt);

Ito_ko = Kko.A3FF;
Ito_ko = nanmean(Ito_ko);
Ito_wt = Kwt.A3;
Ito_wt = nanmean(Ito_wt);

Iss_ko = nanmean(Kko.IssFF);
Iss_wt = nanmean(Kwt.Iss);

tau2_ko = Kko.Tau2FF;
tau2_ko = nanmean(tau2_ko);
tau2_wt = Kwt.Tau2;
tau2_wt = nanmean(tau2_wt);

tau1_ko = Kko.Tau1FF;
tau1_ko = nanmean(tau1_ko);
tau1_wt = Kwt.Tau1;
tau1_wt = nanmean(tau1_wt);

tau3_ko = Kko.Tau3FF;
tau3_ko = nanmean(tau3_ko);
tau3_wt = Kwt.Tau3;
tau3_wt = nanmean(tau3_wt);

t = Ktrace_sub.time;

IKslow1_ko_trace = standard_exp(t, IKslow1_ko, tau2_ko);
IKslow1_ko_trace = transpose(IKslow1_ko_trace);
IKslow1_wt_trace = standard_exp(t, IKslow1_wt, tau2_wt);
IKslow1_wt_trace = transpose(IKslow1_wt_trace);

IKslow2_ko_trace = standard_exp(t, IKslow2_ko, tau1_ko);
IKslow2_ko_trace = transpose(IKslow2_ko_trace);
IKslow2_wt_trace = standard_exp(t, IKslow2_wt, tau1_wt);
IKslow2_wt_trace = transpose(IKslow2_wt_trace);

Ito_ko_trace = standard_exp(t, Ito_ko, tau3_ko);
Ito_ko_trace = transpose(Ito_ko_trace);
Ito_wt_trace = standard_exp(t, Ito_wt, tau3_wt);
Ito_wt_trace = transpose(Ito_wt_trace);


%% draw traces
cap_ko = Kko.CapFF;
cap_ko = nanmean(cap_ko);
cap_wt = Kwt.Cap;
cap_wt = nanmean(cap_wt);

norm_Ktrace_ko = Ktrace_sub.KO ./ cap_ko;
norm_Ktrace_wt = Ktrace_sub.WT ./ cap_wt;

raw_IKslow1_ko_trace = norm_Ktrace_ko - IKslow2_ko_trace - Ito_ko_trace - Iss_ko;
raw_IKslow1_wt_trace = norm_Ktrace_wt - IKslow2_wt_trace - Ito_wt_trace - Iss_wt;

plot(t, raw_IKslow1_ko_trace, 'LineWidth', 2)
hold on
plot(t, raw_IKslow1_wt_trace, 'LineWidth', 2)
hold off
legend('KO', 'WT')


%% investigated fitted parameters - GA IKslow1
% import GA_IKslow1
holding_p = -70; %mV
holding_t = 45; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
P2 = -70; % mV
P2_t = P1_t; % ms

ga_IKslow1_ko = readtable('./results/IKslow1_20200424/IKslow1_KO.csv');
ga_IKslow1_wt = readtable('./results/IKslow1_20200424/IKslow1_WT.csv');
ga_Ito_wt = readtable('./results/GA_Ito_wt_10.xlsx');

for i=1:10
    params = table2array(ga_IKslow1_wt(i, 1:7));
    [t,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,params);
    hold on
    plot(t, A(:,65))
    hold off
end
legend('1','2','3','4','5','6','7','8','9','10')


%% compare traces
IKslow1_ko = mean(table2array(ga_IKslow1_ko(:, 1:7)));
IKslow1_wt = mean(table2array(ga_IKslow1_wt(:, 1:7)));
IKslow2_wt = mean(rst(:, 1:7));
Ito_wt = mean(table2array(ga_Ito_wt(:, 1:6)));

[t1,~,A1,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,IKslow1_wt);
[t2,~,A2,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,IKslow2_wt);
[t3,~,A4,~] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,Ito_wt);

tic
[t,~,A,~] = RasmussonUnparam(holding_p,holding_t,P1,P1_t,P2,P2_t);
toc

figure(1)
plot(t, A(:,65))
hold on
plot(t1, A1(:,65))
plot(t2, A2(:,65))
hold off
legend('Rasmusson', 'KO', 'WT')

[peak, peak_idx] = max(A(:,65));
tau_i = peak*exp(-1);
i_after_peak = A(peak_idx+1:end, 65);
[tau_i_hat, tau_hat] = min(abs(tau_i - i_after_peak));


%% IKsum
