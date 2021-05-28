%% KO - expoenetial fn vs. my model
clc
close all
clear variables

% fitted on 2020-08-28 by SGBA
to_param = [-3.10902175284165,-11.4945930320253,69.4491868862909,37.5847380960043,29.9558997795475,0.170744990408810];
kslow_param = [41.6580470603940,5.47909573395254,31.3763688930612,9.45978973170458,1111.62338230588,0.0227149556759817];
param_ko = [to_param, kslow_param];

% [t, S, A, ~] = AP_param2(param_ko, -91.1, 90);
% plot(t, S(:,1))

holding_p = -70; %mV
holding_t = 0.125*1000; %ms
P1 = 50; %mV
P1_t = 4.5*1000; % ms
Ek = -91.1;

[t1, ~, A1] = Ito(to_param, holding_p, holding_t, P1, P1_t, Ek);
[t2, ~, A2] = IKslow(kslow_param, holding_p, holding_t, P1, P1_t, Ek);

% Ito = A1(:,5);
% IKslow1 = A2(:,5);
% IKslow2 = A3(:,5);

% figure(1)
% plot(t1, Ito)
% hold on
% plot(t2, IKslow1)
% hold off
% axis tight

amps = [24.9, 20.5];
taus = [110.4, 1441.6];

t_holding = 0:0.1:holding_t;
t_P1 = 0:0.1:(P1_t-holding_t);
t_P1_shift = t_P1 + holding_t + 0.1;
t = [t_holding, t_P1_shift];

holding_exp = zeros(1, length(t_holding)); 
Ito_exp = exp_fn(t_P1, amps(1), taus(1));
Ito_exp = [holding_exp, Ito_exp];

IKslow_exp = exp_fn(t_P1, amps(2), taus(2));
IKslow_exp = [holding_exp, IKslow_exp];

% figure(2)
% plot(t, Ito_exp)
% hold on
% plot(t, IKslow_exp)
% hold off
% axis tight

figure(3)
subplot(1,2,1)
plot(t1, A1(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(1,2,2)
plot(t2, A2(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')


%% WT - expoenetial fn vs. my model
% fitted on 2020-08-28 by SGBA
to_param_wt = [36.0323917855634,41.0671435682954,54.3892960635912,53.4110530337884,34.0167375236982,0.224449754332619];
kslow_param_wt = [17.1786204850059,3.55360447475750,32.6483194442754,1.25738530988108,1117.69914528303,0.123107346569954];
param_wt = [to_param_wt, kslow_param_wt];

% AP
% [t, S, A, ~] = AP_param2(param_wt, -91.1, 90);
% plot(t, S(:,1))

% [tt, SS, AA, ~] = Rasmusson_AP(70);
% plot(tt, SS(:,1))

% Kv currents
holding_p = -70; %mV
holding_t = 0.125*1000; %ms
P1 = 50; %mV
P1_t = 4.5*1000; % ms
Ek = -91.1;

[t1_wt, ~, A1_wt] = Ito(to_param_wt, holding_p, holding_t, P1, P1_t, Ek);
[t2_wt, ~, A2_wt] = IKslow(kslow_param_wt, holding_p, holding_t, P1, P1_t, Ek);

Ito_wt = A1_wt(:,5);
IKslow_wt = A2_wt(:,5);

amps_wt = [24.8, 17.1];
taus_wt = [105.2, 1119.6];

t_holding = 0:0.1:holding_t;
t_P1 = 0:0.1:(P1_t-holding_t);
t_P1_shift = t_P1 + holding_t + 0.1;
t = [t_holding, t_P1_shift];
holding_exp = zeros(1, length(t_holding)); 

Ito_exp_wt = exp_fn(t_P1, amps_wt(1), taus_wt(1));
Ito_exp_wt = [holding_exp, Ito_exp_wt];

IKslow_exp_wt = exp_fn(t_P1, amps_wt(2), taus_wt(2));
IKslow_exp_wt = [holding_exp, IKslow_exp_wt];

figure
title('WT')
subplot(1,2,1)
plot(t1_wt, Ito_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
title('I_{to} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(1,2,2)
plot(t2_wt, IKslow_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow1} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')


%% Figure 5 for justifications
figure
subplot(2,2,1)
plot(t1, A1(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
title('I_{to} - KO')
legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,2,2)
plot(t2, A2(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow} - KO')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,2,3)
plot(t1_wt, Ito_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
title('I_{to} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,2,4)
plot(t2_wt, IKslow_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')
