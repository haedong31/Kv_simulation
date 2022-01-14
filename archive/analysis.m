clc
close all
clear variables

holding_p = -70; %mV
holding_t = 450; %ms
P1_t = 25*1000; % ms
Ek = -91.1;

to_param_ko = [-3.10902175284165,-11.4945930320253,68.44919,37.5847380960043,29.9558997795475,0.170744990408810];
kslow1_param_ko = [41.6580470603940,5.47909573395254,31.3763688930612,9.45978973170458,1111.62338230588,0.0227149556759817];
kslow2_param_ko = [41.5408310196983,11.4792079724439,42.0528583047911,8.66747473211152,11263.9824825858,0.0265283219824804];
param_ko = [to_param_ko, kslow1_param_ko, kslow2_param_ko];

to_param_wt = [36.0323917855634,41.0671435682954,54.3892960635912,53.4110530337884,34.0167375236982,0.224449754332619];
kslow1_param_wt = [17.1786204850059,3.55360447475750,32.6483194442754,1.25738530988108,1117.69914528303,0.123107346569954];
kslow2_param_wt = [-23.7191643973443,5.28511873412509,31.7611334578887,14.1815420635500,7231.18383776139,0.0511580416798235];
param_wt = [to_param_wt, kslow1_param_wt, kslow2_param_wt];

%% run
V = -70:10:50;
numv = length(V);

to_tko_ctl = cell(1, numv);
to_Ako_ctl = cell(1, numv);
to_twt_ctl = cell(1, numv);
to_Awt_ctl = cell(1, numv);
to_peaks_ko = zeros(1, numv);
to_peaks_wt = zeros(1, numv);
to_ssa_ko = zeros(1, numv);
to_ssa_wt = zeros(1, numv);

kslow1_tko_ctl = cell(1, numv);
kslow1_Ako_ctl = cell(1, numv);
kslow1_twt_ctl = cell(1, numv);
kslow1_Awt_ctl = cell(1, numv);
kslow1_peaks_ko = zeros(1, numv);
kslow1_peaks_wt = zeros(1, numv);
kslow1_ssa_ko = zeros(1, numv);
kslow1_ssa_wt = zeros(1, numv);

kslow2_tko_ctl = cell(1, numv);
kslow2_Ako_ctl = cell(1, numv);
kslow2_twt_ctl = cell(1, numv);
kslow2_Awt_ctl = cell(1, numv);
kslow2_peaks_ko = zeros(1, numv);
kslow2_peaks_wt = zeros(1, numv);
kslow2_ssa_ko = zeros(1, numv);
kslow2_ssa_wt = zeros(1, numv);
for i=1:numv
    % running simulation
    [to_tko, ~, to_Ako] = Ito(to_param_ko, holding_p, holding_t, V(i), P1_t, Ek);
    [to_twt, ~, to_Awt] = Ito(to_param_wt, holding_p, holding_t, V(i), P1_t, Ek);
    [kslow1_tko, ~, kslow1_Ako] = IKslow(kslow1_param_ko, holding_p, holding_t, V(i), P1_t, Ek);
    [kslow1_twt, ~, kslow1_Awt] = IKslow(kslow1_param_wt, holding_p, holding_t, V(i), P1_t, Ek);
    [kslow2_tko, ~, kslow2_Ako] = IKslow(kslow2_param_ko, holding_p, holding_t, V(i), P1_t, Ek);
    [kslow2_twt, ~, kslow2_Awt] = IKslow(kslow2_param_wt, holding_p, holding_t, V(i), P1_t, Ek);

    % save results
    to_tko_ctl{i} = to_tko;
    to_Ako_ctl{i} = to_Ako;
    to_twt_ctl{i} = to_twt;
    to_Awt_ctl{i} = to_Awt;
    
    kslow1_tko_ctl{i} = kslow1_tko;
    kslow1_Ako_ctl{i} = kslow1_Ako;
    kslow1_twt_ctl{i} = kslow1_twt;
    kslow1_Awt_ctl{i} = kslow1_Awt;

    kslow2_tko_ctl{i} = kslow2_tko;
    kslow2_Ako_ctl{i} = kslow2_Ako;
    kslow2_twt_ctl{i} = kslow2_twt;
    kslow2_Awt_ctl{i} = kslow2_Awt;

    % extract peak information (current densities)
    to_peaks_ko(i) = max(to_Ako(:,5));
    to_peaks_wt(i) = max(to_Awt(:,5));

    kslow1_peaks_ko(i) = max(kslow1_Ako(:,5));
    kslow1_peaks_wt(i) = max(kslow1_Awt(:,5));

    kslow2_peaks_ko(i) = max(kslow2_Ako(:,5));
    kslow2_peaks_wt(i) = max(kslow2_Awt(:,5));

    % extracy SSA information
    to_ssa_ko(i) = to_Ako(end,1) / (to_Ako(end,1) + to_Ako(end,2));
    to_ssa_wt(i) = to_Awt(end,1) / (to_Awt(end,1) + to_Awt(end,2));

    kslow1_ssa_ko(i) = kslow1_Ako(end,1);
    kslow1_ssa_wt(i) = kslow1_Awt(end,1);

    kslow2_ssa_ko(i) = kslow2_Ako(end,1);
    kslow2_ssa_wt(i) = kslow2_Awt(end,1);
end

%% figure; current traces 2 X 3 (group X component) figure
figure(1)
set(gca,'Box','on');

subplot(2,3,1)
title('I_{to} - KO')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
for i=1:numv
    hold on
        plot(to_tko_ctl{i}, to_Ako_ctl{i}(:,5))
    hold off
end
axis tight

subplot(2,3,2)
title('I_{Kslow1} - KO')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
for i=1:numv
    hold on
        plot(kslow1_tko_ctl{i}, kslow1_Ako_ctl{i}(:,5))
    hold off
end
axis tight

subplot(2,3,3)
title('I_{Kslow2} - KO')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
for i=1:numv
    hold on
        plot(kslow2_tko_ctl{i}, kslow2_Ako_ctl{i}(:,5))
    hold off
end
axis tight

subplot(2,3,4)
title('I_{to} - WT')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
for i=1:numv
    hold on
        plot(to_twt_ctl{i}, to_Awt_ctl{i}(:,5))
    hold off
end
axis tight

subplot(2,3,5)
title('I_{Kslow1} - WT')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
for i=1:numv
    hold on
        plot(kslow1_twt_ctl{i}, kslow1_Awt_ctl{i}(:,5))
    hold off
end
axis tight

subplot(2,3,6)
title('I_{Kslow2} - WT')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
for i=1:numv
    hold on
        plot(kslow2_twt_ctl{i}, kslow2_Awt_ctl{i}(:,5))
    hold off
end
axis tight

%% figure; current densities 1 X 3 horizontal figure
figure(2)
subplot(1,3,1)
xlabel('Voltage (mV)')
ylabel('A_{to} (pA/pF)')
hold on
    plot(V, to_peaks_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, to_peaks_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
axis tight
grid on
legend('MGAT1KO', 'WT')
xticks(V)

subplot(1,3,2)
xlabel('Voltage (mV)')
ylabel('A_{Kslow1} (pA/pF)')
hold on
    plot(V, kslow1_peaks_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow1_peaks_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
axis tight
grid on
% legend('MGAT1KO', 'WT')
xticks(V)

subplot(1,3,3)
xlabel('Voltage (mV)')
ylabel('A_{Kslow2} (pA/pF)')
hold on
    plot(V, kslow2_peaks_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow2_peaks_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
axis tight
grid on
% legend('MGAT1KO', 'WT')
xticks(V)

%% figure; SSA 1 X 3 horizontal figure
figure(3)
subplot(1,3,1)
xlabel('Voltage (mV)')
ylabel('SSA')
hold on
    plot(V, to_ssa_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, to_ssa_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
axis tight
grid on
legend('MGAT1KO', 'WT')
xticks(V)
ylim([0, 1])

subplot(1,3,2)
xlabel('Voltage (mV)')
ylabel('SSA')
hold on
    plot(V, kslow1_ssa_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow1_ssa_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
axis tight
grid on
% legend('MGAT1KO', 'WT')
xticks(V)
ylim([0, 1])

subplot(1,3,3)
xlabel('Voltage (mV)')
ylabel('SSA')
hold on
    plot(V, kslow2_ssa_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow2_ssa_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
axis tight
grid on
% legend('MGAT1KO', 'WT')
xticks(V)
ylim([0, 1])
