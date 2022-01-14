clc
close all
clear variables


%% KO
K_ko = readtable("./potassium-KO.xlsx");
IKslow2_ko = K_ko.A1FF;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [0.00000481333, 26.5, 26.5, 0.0000953333, 26.5]
num_vars = 5;
lower_bd = [-0.00000481333, -26.5, -26.5, -0.0000953333, -26.5];
upper_bd = [1.443999e-05, 79.5, 79.5, 0.0002859999, 79.5];
fit_fn = @(X) IKslow2_fitness(X,IKslow2_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 10);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow2(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow2 = A(:,64);
    
    % save the results
    file_path = sprintf('./GA_IKslow2_ko_%i.mat', i);
    rst(i,:) = [X, fval, fitted_IKslow2(end)];
    save(file_path, 'rst');
    disp(rst)
end


%% WT
K_wt = readtable("./potassium-WT.xlsx");
IKslow2_wt = K_wt.A1;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [0.00000481333, 26.5, 26.5, 0.0000953333, 26.5]
num_vars = 5;
lower_bd = [-0.00000481333, -26.5, -26.5, -0.0000953333, -26.5];
upper_bd = [1.443999e-05, 79.5, 79.5, 0.0002859999, 79.5];
fit_fn = @(X) IKslow2_fitness(X,IKslow2_wt,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 10);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow2(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow2 = A(:,64);
    
    % save the results
    file_path = sprintf('./GA_IKslow2_wt_%i.mat', i);
    rst(i,:) = [X, fval, fitted_IKslow2(end)];
    save(file_path, 'rst');
    disp(rst)
end
