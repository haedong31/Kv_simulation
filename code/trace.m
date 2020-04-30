%% manipulate the average raw trace
clc
close all
clear variables

Ktrace = readtable('k_trace.csv');
Ktrace_sub = Ktrace(9861:500001,:);
Ktrace_sub.time = Ktrace_sub.time - 493;

% plot(Ktrace_sub.time, Ktrace_sub.KO)


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


%%
plot(t, Ito_ko_trace, 'LineWidth', 2)
hold on
plot(t, Ito_wt_trace, 'LineWidth', 2)
hold off
legend('KO', 'WT')
