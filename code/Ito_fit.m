clc
close all
clear variables


%% Ito - KO
% import data
K_ko = readtable("./potassium-KO.xlsx");
Ito_ko = K_ko.A3FF;

% voltage clamp protocol parameters
holding_p = -70; %mV
holding_t = 4.5; %msec
P1 = 50; %mV
P1_t = 29.5; % msec
P2 = 50; % mV
P2_t = 29.5; % msec

% GA parameters
% X = [30,30,13.5,33.5,33.5,33.5];
num_vars = 6;
lower_bd = [0,0,0,0,0,0];
upper_bd = [60,60,60,60,60,60];
fit_fn = @(X) Ito_fitness(X,Ito_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 8);
for i=1:10
    tic
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value (steady-state current)
    fitted_Ito = A(:,61);
    
    % save the results
    file_path = sprintf('./GA_Ito_ko_%i.mat', i);
    rst(i,:) = [X, fval, max(fitted_Ito)];
    save(file_path, 'rst');
    disp(rst)
    toc
end


%% Ito - WT
K_wt = readtable("./potassium-WT.xlsx");
Ito_wt = K_wt.A3;

% voltage clamp protocol parameters
holding_p = -70; %mV
holding_t = 4.5; %msec
P1 = 50; %mV
P1_t = 29.5; % msec
P2 = 50; % mV
P2_t = 29.5; % msec

% GA parameters
num_vars = 6;
lower_bd = [0,0,0,0,0,0];
upper_bd = [60,60,60,60,60,60];
fit_fn = @(X) Ito_fitness(X,Ito_wt,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 8);
for i=1:10
    tic
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value (steady-state current)
    fitted_Ito = A(:,61);
    
    % save the results
    file_path = sprintf('./GA_Ito_wt_%i.mat', i);
    rst(i,:) = [X, fval, max(fitted_Ito)];
    save(file_path, 'rst');
    disp(rst)
    toc
end
