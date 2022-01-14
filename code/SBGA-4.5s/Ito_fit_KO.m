%% KO
clc
close all
clear variables
format shortG

% data to fit model
y = [18.39, 111.86];
tol = [0.1, 1.0];

% model fitting
num_iters = 30;
best_chroms = zeros(num_iters, 6);
best_amps = zeros(1, num_iters);
best_taus = zeros(1, num_iters);
for i = 1:num_iters
    fprintf('### Iteration %i \n', i)
    [amps, taus, ~, chroms] = Ito_SBGA(y, tol, 30, 6, 4);
    best_chroms(i, :) = chroms(end, :);
    best_amps(i) = amps(end);
    best_taus(i) = taus(end);
    clc
end

save('Ito_KO.mat', 'best_chroms', 'best_amps', 'best_taus')
