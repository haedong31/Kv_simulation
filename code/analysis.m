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
% import IKslow1_ko_ofWhatYouWant
% import IKslow1_wt_ofWhatYouWant
ga_ko = mean(table2array(IKslow1KO));
ga_wt = mean(table2array(IKslow1WT));
rt_ko = [-68.0570165036410	2.58308797860332	50.8876752241103	1434.81281801083	18.0013636837798	0.849759529744222	35.2985494078469];
rt_wt = [24.1679224667629	5.56468271665088	66.0283391841861	2316.27619578639	62.8289892211523	1.46390048491702	136.674046649574];

[t1, ~, A1, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,ga_ko);
[t2, ~, A2, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,ga_wt);
[t3, ~, A3, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,rt_ko);
[t4, ~, A4, ~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,rt_wt);

plot(t1, A1(:,65))
hold on
plot(t2, A2(:,65))
plot(t3, A3(:,65))
plot(t4, A4(:,65))
hold off
legend('GA KO', 'GA WT', 'RT KO', 'RT WT')


%% simple_RF analysis
plot(deltas)
xlabel('Iteration')
ylabel('Model Discrepancy')
title('IKslow1 KO')

[d, d_idx] = min(delta);
p = params(d_idx,:);
