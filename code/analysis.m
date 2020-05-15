clc
close all
clear variables


%% K+ protocol
holding_p = -70; %mV
holding_t = 4.5*10; %ms
P1 = 50; %mV
P1_t = 29.5*10; % ms
P2 = -70; % mV
P2_t = P1_t; % ms


%% Ito - run simulation
load('./results/GA_Ito_ko_10.mat')
ga_ito_ko = mean(rst, 1);
load('./results/GA_Ito_wt_10.mat')
ga_ito_wt = mean(rst, 1);
Ito_trace = readtable('./Ito_trace.csv');

[t1, ~, A1, ~] = Ito(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_ito_ko(1:6));
[t2, ~, A2, ~] = Ito(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_ito_wt(1:6));

plot(Ito_trace.time, Ito_trace.KO)
hold on
plot(Ito_trace.time, Ito_trace.WT)
hold off
legend('KO', 'WT')
title('Average Ito Trace')


%% draw traces
figure(1)
t = 0:0.00326:296;
t(90362:end) = [];
yyaxis left
plot(t1, A1(:,61)*257.9375)
yyaxis right
plot(t, Ito_trace.KO(9641:100001))
title('KO')
legend('Simulated', 'Real Data')

% Wt
figure(2)
yyaxis left
plot(t2, A2(:,61)*257.9375)
yyaxis right
plot(t, Ito_trace.WT(9641:100001))
title('WT')
legend('Simulated', 'Real Data')

% simulated results
figure(3)
plot(t1, A1(:,61)*257.9375)
hold on
plot(t2, A2(:,61)*257.9375)
hold off
legend('KO', 'WT')
title('Simulated Ito Traces')


%% Iss - Compare the simulated and raw traces
load('./results/GA_Iss_ko_10.mat')
ga_iss_ko = mean(rst, 1);
load('./results/GA_Iss_wt_10.mat')
ga_iss_wt = mean(rst, 1);

[t3, ~, A3, ~] = Iss(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_iss_ko(1:4));
[t4, ~, A4, ~] = Iss(holding_p, holding_t, P1, P1_t, P2, P2_t, ga_iss_wt(1:4));


%% visualization
% import ga_ikslow1_ko
% import ga_ikslow1_wt

X_ko = mean(IKslow1KO);
X_wt = mean(IKslow1WT);

[t, S, A, ~] = Rasmusson_AP(1000);
[t1, S1, A1, ~] = IKslow1_AP(1000, X_ko(1:7));
[t2, S2, A2, ~] = IKslow1_AP(1000, IKslow1WT(2,1:7));

plot(t1, S1(:,1))
hold on
plot(t2, S2(:,1))
hold off
legend('KO', 'WT')


%% simple_RF analysis
plot(deltas)
xlabel('Iteration')
ylabel('Model Discrepancy')

[d, d_idx] = min(delta);
p = params(d_idx,:);
