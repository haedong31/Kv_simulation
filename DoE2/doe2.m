clc
close all
clear variables


%% import data
cap_ko = 254.3;

raw_trc = readtable('k_trace.csv');
raw_trc = raw_trc(1:500001,:);

raw_ko = raw_trc.KO;
raw_ko = raw_ko ./ cap_ko;
raw_ko(raw_ko < 0) = 0;
raw_time = raw_trc.time;

% chop off early phase
[~, peak_idx] = max(raw_ko);
raw_ko_rd = raw_ko(peak_idx:end);
raw_time_rd = raw_time(peak_idx:end);
raw_time_rd = raw_time_rd - raw_time_rd(1);

raw_trc_rd = table(raw_time_rd, raw_ko_rd);
raw_trc_rd.Properties.VariableNames = {'time', 'current'};
ds_trc = downsample(raw_trc_rd, 400);


%% 
V = -80:1:60;
x9 = -81.46;
x10 = -5.98;
x11 = -34.31;
x12 = -4.31;
irs1 = 0.21 ./ (1+exp(-(V-x9)./x10));
irs2 = 0.79 ./ (1+exp(-(V-x11)./x12));
plot(V, irs1, 'LineWidth',2)
hold on
plot(V, irs2, 'LineWidth',2)
hold off
legend('irs1','irs2')
% irs = 0.21 ./ (1+exp(-(V-x9)./x10)) + 0.79 ./ (1+exp(-(V-x11)./x12));


