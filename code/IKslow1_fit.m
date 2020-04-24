clc
close all
clear variables


%% KO
K_ko = readtable("./potassium-KO.xlsx");
IKslow1_ko = K_ko.A2FF;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [22.5000, 2.05800, 45.2000, 1200.00, 45.2000, 0.493000, 170.000]
num_vars = 7;
lower_bd = [-22.5, -2.058, -45.2, -1200.0, -45.0, -0.5, -100.0];
upper_bd = [67.5, 6.174, 135.6, 3600.0, 90.0, 1.5, 400.0];
fit_fn = @(X) IKslow1_fitness(X,IKslow1_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 9);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow1 = A(:,65);
    
    % save the results
    file_path = sprintf('./GA_IKslow1_ko_%i.mat', i);
    rst(i,:) = [X, fval, max(fitted_IKslow1)];
    save(file_path, 'rst');
    disp(rst)
end


%% WT
K_wt = readtable("./potassium-WT.xlsx");
IKslow1_wt = K_wt.A2;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [22.5000, 2.05800, 45.2000, 1200.00, 45.2000, 0.493000, 170.000]
num_vars = 7;
lower_bd = [-22.5, -2.058, -45.2, -1200.0, -45.0, -0.5, -100.0];
upper_bd = [67.5, 6.174, 135.6, 3600.0, 90.0, 1.5, 400.0];
fit_fn = @(X) IKslow1_fitness(X,IKslow1_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 9);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow1 = A(:,65);
    
    % save the results
    file_path = sprintf('./GA_IKslow1_wt_%i.mat', i);
    rst(i,:) = [X, fval, max(fitted_IKslow1)];
    save(file_path, 'rst');
    disp(rst)
end
