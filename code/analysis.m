clc
close all
clear variables


%% current densities and taus
holding_p = -70; %mV
holding_t = 450; %ms
P1 = -60:10:50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

init_to_ko = [-13.5655  128.4098  321.7877  127.2189   58.4796];
init_Kslow1_ko = [-0.0613    0.0097    0.2070    0.0128    1.1628];
init_Kslow1_ko = init_Kslow1_ko*1000;
init_Kslow2_ko = [-0.0717    0.0123    0.0245    0.0399    8.6985];
init_Kslow2_ko = init_Kslow2_ko*1000;
init_param_ko = [init_to_ko init_Kslow1_ko init_Kslow2_ko];

init_to_wt = [0.661280350042337  -0.706849509040702   1.007979522932042   0.813067763455725   0.472446214481949];
init_to_wt = init_to_wt*100;
init_Kslow1_wt = [0.015413685966188   0.008810185337584   0.076048262113858   0.018025278517100   1.150596392908760];
init_Kslow1_wt = init_Kslow1_wt*1000;
init_Kslow2_wt = [-0.055151382655153   0.006381105415206  -0.025994281704011   0.018749874829633   4.424133917575417];
init_Kslow2_wt = init_Kslow2_wt*1000;
init_param_wt = [init_to_wt init_Kslow1_wt init_Kslow2_wt];

to_peaks_ko = zeros(1, length(P1));
to_peaks_wt = zeros(1, length(P1));
Kslow1_peaks_ko = zeros(1, length(P1));
Kslow1_peaks_wt = zeros(1, length(P1));
Kslow2_peaks_ko = zeros(1, length(P1));
Kslow2_peaks_wt = zeros(1, length(P1));

to_taus_ko = zeros(1, length(P1));
to_taus_wt = zeros(1, length(P1));
Kslow1_taus_ko = zeros(1, length(P1));
Kslow1_taus_wt = zeros(1, length(P1));
Kslow2_taus_ko = zeros(1, length(P1));
Kslow2_taus_wt = zeros(1, length(P1));

for i=1:length(P1)
    [t1, ~, A1, ~] = IKsum(init_param_ko, holding_p, holding_t, P1(i), P1_t, Ek);
    [t2, ~, A2, ~] = IKsum(init_param_wt, holding_p, holding_t, P1(i), P1_t, Ek);
    
    Ito_trc_ko = A1(:,5);
    Ito_trc_wt = A2(:,5);
    IKslow1_trc_ko = A1(:,10);
    IKslow1_trc_wt = A2(:,10);
    IKslow2_trc_ko = A1(:,15);
    IKslow2_trc_wt = A2(:,15);

    [to_peak_ko, to_peak_ko_idx] = max(Ito_trc_ko);
    [to_peak_wt, to_peak_wt_idx] = max(Ito_trc_wt);
    [Kslow1_peak_ko, Kslow1_peak_ko_idx] = max(IKslow1_trc_ko);
    [Kslow1_peak_wt, Kslow1_peak_wt_idx] = max(IKslow1_trc_wt);
    [Kslow2_peak_ko, Kslow2_peak_ko_idx] = max(IKslow2_trc_ko);
    [Kslow2_peak_wt, Kslow2_peak_wt_idx] = max(IKslow2_trc_wt);

    to_peaks_ko(i) = to_peak_ko;
    to_peaks_wt(i) = to_peak_wt;
    Kslow1_peaks_ko(i) = Kslow1_peak_ko;
    Kslow1_peaks_wt(i) = Kslow1_peak_wt;
    Kslow2_peaks_ko(i) = Kslow2_peak_ko;
    Kslow2_peaks_wt(i) = Kslow2_peak_wt;

    [~, to_tau_ko_idx] = min(abs(to_peak_ko*exp(-1) - Ito_trc_ko(to_peak_ko_idx:end)));
    [~, to_tau_wt_idx] = min(abs(to_peak_wt*exp(-1) - Ito_trc_wt(to_peak_wt_idx:end)));
    [~, Kslow1_tau_ko_idx] = min(abs(Kslow1_peak_ko*exp(-1) - IKslow1_trc_ko(Kslow1_peak_ko_idx:end)));
    [~, Kslow1_tau_wt_idx] = min(abs(Kslow1_peak_wt*exp(-1) - IKslow1_trc_wt(Kslow1_peak_wt_idx:end)));
    [~, Kslow2_tau_ko_idx] = min(abs(Kslow2_peak_ko*exp(-1) - IKslow2_trc_ko(Kslow2_peak_ko_idx:end)));
    [~, Kslow2_tau_wt_idx] = min(abs(Kslow2_peak_wt*exp(-1) - IKslow2_trc_wt(Kslow2_peak_wt_idx:end)));

    to_taus_ko = t1(to_tau_ko_idx);
    to_taus_wt = t2(to_tau_wt_idx);
    Kslow1_taus_ko = t1(Kslow1_tau_ko_idx);
    Kslow1_taus_wt = t2(Kslow1_tau_wt_idx);
    Kslow2_taus_ko = t1(Kslow2_tau_ko_idx);
    Kslow2_taus_wt = t2(Kslow2_tau_wt_idx);
end

% to_max_peak_ko = max(to_peaks_ko);
% to_max_peak_wt = max(to_peaks_wt);
% to_peaks_ko = to_peaks_ko ./ to_max_peak_ko;
% to_peaks_wt = to_peaks_wt ./ to_max_peak_wt;
% 
% Kslow1_max_peak_ko = max(Kslow1_peaks_ko);
% Kslow1_max_peak_wt = max(Kslow1_peaks_wt);
% Kslow1_peaks_ko = Kslow1_peaks_ko ./ Kslow1_max_peak_ko;
% Kslow1_peaks_wt = Kslow1_peaks_wt ./ Kslow1_max_peak_wt;
% 
% Kslow2_max_peak_ko = max(Kslow2_peaks_ko);
% Kslow2_max_peak_wt = max(Kslow2_peaks_wt);
% Kslow2_peaks_ko = Kslow2_peaks_ko ./ Kslow2_max_peak_ko;
% Kslow2_peaks_wt = Kslow2_peaks_wt ./ Kslow2_max_peak_wt;

figure(1)
plot(P1, to_peaks_ko, '-o', 'Color','red')
hold on
plot(P1, to_peaks_wt, '-*', 'Color','blue')
hold off
axis tight
legend('KO','WT')
ylabel('pA/pF')
xlabel('mV')
title('Current Density Ito')

figure(2)
plot(P1, Kslow1_peaks_ko, '-o', 'Color','red')
hold on
plot(P1, Kslow1_peaks_wt, '-*', 'Color','blue')
hold off
axis tight
legend('KO','WT')
ylabel('pA/pF')
xlabel('mV')
title('Current Density IKslow1')

figure(3)
plot(P1, Kslow2_peaks_ko, '-o', 'Color','red')
hold on
plot(P1, Kslow2_peaks_wt, '-*', 'Color','blue')
hold off
axis tight
legend('KO','WT')
ylabel('pA/pF')
xlabel('mV')
title('Current Density IKslow2')

figure(4)
plot(P1, to_taus_ko, '-o', 'Color','red')
hold on
plot(P1, to_taus_wt, '-*', 'Color','blue')
hold off
legend('KO','WT')
ylabel('ms')
xlabel('mV')
title('Time constant Ito')

figure(5)
plot(P1, Kslow1_taus_ko, '-o', 'Color','red')
hold on
plot(P1, Kslow1_taus_wt, '-*', 'Color','blue')
hold off
legend('KO','WT')
ylabel('ms')
xlabel('mV')
title('Time constant IKslow1')

figure(6)
plot(P1, Kslow2_taus_ko, '-o', 'Color','red')
hold on
plot(P1, Kslow2_taus_wt, '-*', 'Color','blue')
hold off
legend('KO','WT')
ylabel('ms')
xlabel('mV')
title('Time constant IKslow2')


% %% Ito - run simulation
% load('./results/GA_Ito_ko_10.mat')
% ga_ito_ko = mean(rst, 1);
% load('./results/GA_Ito_wt_10.mat')
% ga_ito_wt = mean(rst, 1);
% Ito_trace = readtable('./Ito_trace.csv');
% 
% [t1, ~, A1, ~] = Ito(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_ito_ko(1:6));
% [t2, ~, A2, ~] = Ito(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_ito_wt(1:6));
% 
% plot(Ito_trace.time, Ito_trace.KO)
% hold on
% plot(Ito_trace.time, Ito_trace.WT)
% hold off
% legend('KO', 'WT')
% title('Average Ito Trace')
% 
% 
% %% draw traces
% figure(1)
% t = 0:0.00326:296;
% t(90362:end) = [];
% yyaxis left
% plot(t1, A1(:,61)*257.9375)
% yyaxis right
% plot(t, Ito_trace.KO(9641:100001))
% title('KO')
% legend('Simulated', 'Real Data')
% 
% % Wt
% figure(2)
% yyaxis left
% plot(t2, A2(:,61)*257.9375)
% yyaxis right
% plot(t, Ito_trace.WT(9641:100001))
% title('WT')
% legend('Simulated', 'Real Data')
% 
% % simulated results
% figure(3)
% plot(t1, A1(:,61)*257.9375)
% hold on
% plot(t2, A2(:,61)*257.9375)
% hold off
% legend('KO', 'WT')
% title('Simulated Ito Traces')
% 
% 
% %% Iss - Compare the simulated and raw traces
% load('./results/GA_Iss_ko_10.mat')
% ga_iss_ko = mean(rst, 1);
% load('./results/GA_Iss_wt_10.mat')
% ga_iss_wt = mean(rst, 1);
% 
% [t3, ~, A3, ~] = Iss(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_iss_ko(1:4));
% [t4, ~, A4, ~] = Iss(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_iss_wt(1:4));
% 
% 
% %% visualization
% % import IKslow1_ko_ofWhatYouWant
% % import IKslow1_wt_ofWhatYouWant
% ga_ko = mean(table2array(IKslow1KO));
% ga_wt = mean(table2array(IKslow1WT));
% rt_ko = [-68.0570165036410	2.58308797860332	50.8876752241103	1434.81281801083	18.0013636837798	0.849759529744222	35.2985494078469];
% rt_wt = [-45.4668111833352	5.66962515497214	94.0192673328295	-180.129531987330	18.0619040526542	-0.256041915755219	456.849512275044];
% 
% [t1, ~, A1, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,ga_ko);
% [t2, ~, A2, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,ga_wt);
% [t3, ~, A3, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,rt_ko);
% [t4, ~, A4, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,rt_wt);
% 
% plot(t1, A1(:,65))
% hold on
% plot(t2, A2(:,65))
% plot(t3, A3(:,65))
% plot(t4, A4(:,65))
% hold off
% legend('GA KO', 'GA WT', 'RT KO', 'RT WT')
% 
% 
% %% simple_RF analysis
% plot(min_deltas)
% xlabel('Iteration')
% ylabel('Model Discrepancy')
% title('IKslow1 WT')
% 
% [d, d_idx] = min(delta);
% p = params(d_idx,:);
% 
% [t,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,p);
% plot(t, A(:,65))
