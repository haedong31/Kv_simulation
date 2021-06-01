%% WT
clc
close all
clear variables
format shortG

% data to fit model
y = [4.70, 1874.92];
tol = [0.1, 5.0];

% model fitting
num_iters = 30;
best_chroms = zeros(num_iters, 6);
best_amps = zeros(1, num_iters);
best_taus = zeros(1, num_iters);
for i = 17:num_iters
    fprintf('### Iteration %i \n', i)
    [amps, taus, ~, chroms] = IKslow_SBGA(y, tol, 30, 6, 4);
    best_chroms(i, :) = chroms(end, :);
    best_amps(i) = amps(end);
    best_taus(i) = taus(end);
    clc
end

save('IKslow_KO.mat', 'best_chroms', 'best_amps', 'best_taus')
