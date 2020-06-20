clc
close all
clear variables
format long


% holding_p = -70; %mV
% holding_t = 450; %ms
% P1 = 50; %mV
% P1_t = 25*1000; % ms
% Ek = -91.1;
% X0 = [22.5, 7.7, 45.2, 5.7, 1200.0, 0.16];
% [t, ~, A] = IKslow(X0, holding_p, holding_t, P1, P1_t, Ek);
% plot(t, A(:,5))

df = readtable("./potassium-KO.xlsx");
num_obs = 34;
amp = df.A2;
tau = df.Tau2;

best_amp_container = zeros(1, num_obs);
best_tau_container = zeros(1, num_obs);
num_iters_container = zeros(1, num_obs);
best_chrom_container = zeros(num_obs, 5);
for i=1:num_obs
    if i == 33
        continue
    end

    fprintf('### Iter %i/%i \n', i, num_obs)
    
    y = [amp(i), tau(i)];
    [best_amps, best_taus, best_gens, best_chroms] = IKslow1_AGA(5, y, 30, 6, 4);
    best_amp_container(i) = best_amps(end);
    best_tau_container(i) = best_taus(end);
    num_iters_container(i) = best_gens(end);
    best_chrom_container(i,:) = best_chroms(end,:);
end

% [t1, ~, A1] = IKslow(param_ko(1,:), holding_p, holding_t, P1, P1_t, Ek);
% IKslow1 = A1(:,5);
% [t2, ~, A2] = IKslow(param_ko(2,:), holding_p, holding_t, P1, P1_t, Ek);
% IKslow2 = A2(:,5);
% 
% plot(t1, IKslow1)
% hold on
% plot(t2, IKslow2)
% hold off
% legend('IKslow1','IKslow2')
