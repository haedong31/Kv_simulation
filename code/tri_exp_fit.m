clc
close all
clear variables


%% import data
cap_ko = 254.3;
N0 = 30;
N1 = 6;
N2 = 4;
num_var = 7;

raw_trc = readtable('k_trace.csv');
raw_trc = raw_trc(1:500001,:);

raw_ko = raw_trc.KO;
raw_ko = raw_ko ./ cap_ko;
raw_time = raw_trc.time;

% chop off early phase
[peak, peak_idx] = max(raw_ko);
raw_ko_rd = raw_ko(peak_idx:end);
raw_ko_rd(raw_ko_rd < 0) = 0;
raw_time_rd = raw_time(peak_idx:end);
raw_time_rd = raw_time_rd - raw_time_rd(1);

raw_trc_rd = table(raw_time_rd, raw_ko_rd);
raw_trc_rd.Properties.VariableNames = {'time', 'current'};


%% initialization
% Iss, A1, A2, A3, Tau1, Tau2, Tau3
% A1 - IKslow2; A2 - IKslow1; A3 - Ito
low = [2.08, 1.81, 0.43, 2.13, 4476.04, 262.21,	84.47];
high = [13.50, 6.35, 12.52,	60.99, 57941.35, 2000.97, 203.91];

init_gen = zeros(N0, num_var);
for j=1:num_var
    unif = makedist('Uniform', 'lower',low(j), 'upper',high(j));
    init_gen(:,j) = random(unif, N0, 1);
end


%% evaluation
best_fits = [];

% sampling for evaluation
ds_trc = downsample(raw_trc_rd, 400);
t = ds_trc.time;
fits = zeros(1, N0);
for i=1:N0
    param = init_gen(i,:);
    Iss = param(1);
    IKslow2 = exp_fn(t, param(2), param(5));
    IKslow1 = exp_fn(t, param(3), param(6));
    Ito = exp_fn(t, param(4), param(7));
    
    IKsum = Iss + IKslow2 + IKslow1 + Ito;
    fits(i) = sum((IKsum - ds_trc.current).^2);
end
best_fits = [best_fits, min(fits)];

new_gen = zeros(N0, num_var);
[~, super_idx] = mink(fits, N1);
elites = init_gen(super_idx,:);
new_gen(1:N1,:) = elites;

cnt = 1;
sigs = std(elites);
for i=1:N1
    elite = elites(i);
    for j=1:N2
        offspring = elite + normrnd(0,sigs);
        offspring = abs(offspring);
        new_gen((N1+cnt),:) = offspring;
        cnt = cnt + 1;
    end
end

fits = zeros(1, N0);
for i=1:N0
    param = new_gen(i,:);
    Iss = param(1);
    IKslow2 = exp_fn(t, param(2), param(5));
    IKslow1 = exp_fn(t, param(3), param(6));
    Ito = exp_fn(t, param(4), param(7));
    
    IKsum = Iss + IKslow2 + IKslow1 + Ito;
    fits(i) = sum((IKsum - ds_trc.current).^2);
end
best_fits = [best_fits, min(fits)];

upper_cnt = 1;
for k=1:10000
    fprintf('Generation %i \n', upper_cnt);
    new_gen = zeros(N0, num_var);
    [~, super_idx] = mink(fits, N1);
    elites = init_gen(super_idx,:);
    new_gen(1:N1,:) = elites;

    cnt = 1;
    sigs = std(elites);
    for i=1:N1
        elite = elites(i);
        for j=1:N2
            offspring = elite + normrnd(0,sigs);
            offspring = abs(offspring);
            new_gen((N1+cnt),:) = offspring;
            cnt = cnt + 1;
        end
    end

    fits = zeros(1, N0);
    for i=1:N0
        param = new_gen(i,:);
        Iss = param(1);
        IKslow2 = exp_fn(t, param(2), param(5));
        IKslow1 = exp_fn(t, param(3), param(6));
        Ito = exp_fn(t, param(4), param(7));

        IKsum = Iss + IKslow2 + IKslow1 + Ito;
        fits(i) = sum((IKsum - ds_trc.current).^2);
    end
    best_fits = [best_fits, min(fits)];
    upper_cnt = upper_cnt + 1;
end
