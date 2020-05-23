clc
close all
clear variables
warning('off')


%% data and modeling parameters
% K_ko = readtable('./potassium-KO.xlsx');
% IKslow1 = mean(K_ko.A2FF);
K_wt = readtable('./potassium-WT.xlsx');
IKslow1 = mean(K_wt.A2);

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; 
P1 = 50; % mV
P1_t = 200;
P2 = 50; % mV
P2_t = P1_t;

% original values in the Rasmusson
X = [22.5000, 2.05800, 45.2000, 1200.00, 45.2000, 0.493000, 170.000];
num_params = length(X);
lower_bd = zeros(1, num_params);
upper_bd = zeros(1, num_params);
for i=1:num_params
    lower_bd(i) = X(i) - X(i)*2;
    upper_bd(i) = X(i) + X(i)*2;
end

init_size = 100;
sample_size = 100;
var_rdc_rate = 0.5;
obs_rdc_rate = 0.5;
num_trees = 30;
num_imp_vars = floor(num_params*var_rdc_rate);
num_top_rows = floor(init_size*obs_rdc_rate);
tree_tem = templateTree('NumVariablesToSample','all', 'PredictorSelection','interaction-curvature');
tol = 0.01;


%% initial parameter set & evaluation
% sampling
init_params = zeros(init_size, num_params);
for i=1:num_params
    unif = makedist('Uniform', 'lower',lower_bd(i), 'upper',upper_bd(i));
    init_params(:,i) = random(unif, init_size, 1);
end

% evaluation
IKslow1_hat = zeros(init_size, 1);
for i=1:init_size
    fprintf('%i \n', i)
    [~,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,init_params(i,:));
    IKslow1_hat(i) = max(A(:,65));
end

delta = abs(IKslow1 - IKslow1_hat);
init_data = [init_params, delta];
init_data = array2table(init_data);

% mdl = TreeBagger(30, init_data, 'init_data8', 'Method','regression', 'OOBPredictorImportance','on');
mdl = fitrensemble(init_data, 'init_data8', 'Method','Bag', 'NumLearningCycles',num_trees, 'Learners',tree_tem);
impOOB = oobPermutedPredictorImportance(mdl);

[imp_params, imp_params_idx] = maxk(impOOB, num_imp_vars);
delta = table2array(init_data(:,8));
[imp_rows, imp_rows_idx] = mink(delta, num_top_rows);
nonimp_params_idx = 1:1:num_params;
nonimp_params_idx = setdiff(nonimp_params_idx, imp_params_idx);

% new parameter set
sampled_data = init_data(imp_rows_idx, imp_params_idx);
mu = mean(table2array(sampled_data));
sigma = cov(table2array(sampled_data));
params = zeros(sample_size, num_params);
params(:,imp_params_idx) = mvnrnd(mu, sigma, sample_size);
for i=1:length(nonimp_params_idx)
    idx = nonimp_params_idx(i);
    unif = makedist('Uniform', 'lower',lower_bd(idx), 'upper',upper_bd(idx));
    params(:,idx) = random(unif, sample_size, 1);
end


%% run learning
err_idx = [];
min_deltas = [];
k = 1;
tic
while 1
    tic
    
    fprintf('###Iter %i \n', k);

    % evaluate the proposed parameters with the Rasmusson
    IKslow1_hat = zeros(sample_size, 1);
    for i=1:sample_size
        [~,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,params(i,:));
        IKslow1_hat(i) = max(A(:,65));
    end

    % fit random forest
    delta = abs(IKslow1 - IKslow1_hat);
    new_data = [params, delta];
    new_data = array2table(new_data);

    min_delta = min(delta);
    fprintf('Min delta: %6.4f \n', min_delta)
    min_deltas = [min_deltas; min_delta];
    
    if min_delta <= tol
        break
    end

    % OOB permuted variable importance
    mdl = fitrensemble(new_data, 'new_data8', 'Method','Bag', 'NumLearningCycles',num_trees, 'Learners',tree_tem);
    try
        impOOB = oobPermutedPredictorImportance(mdl);
    catch ME
        err_idx = [err_idx; k];
        fprintf('ERROR at %i: %s \n', k, ME.message)
    end

    % variable selection
    [imp_params, imp_params_idx] = maxk(impOOB, num_imp_vars);
    delta = table2array(new_data(:,8));
    [imp_rows, imp_rows_idx] = mink(delta, num_top_rows);
    nonimp_params_idx = 1:1:num_params;
    nonimp_params_idx = setdiff(nonimp_params_idx, imp_params_idx);

    % mean and covariance for significant parameters
    sampled_data = new_data(imp_rows_idx, imp_params_idx);
    mu = mean(table2array(sampled_data));
    sigma = cov(table2array(sampled_data));
    
    params(:,imp_params_idx) = mvnrnd(mu, sigma, sample_size);
    for i=1:length(nonimp_params_idx)
        idx = nonimp_params_idx(i);
        unif = makedist('Uniform', 'lower',lower_bd(idx), 'upper',upper_bd(idx));
        params(:,idx) = random(unif, sample_size, 1);
    end
    k = k+1;

    toc
end
toc
