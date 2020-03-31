clc
close all
clear variables


%% Iss - KO
% import data
K_ko = readtable("./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx");
Iss_ko = K_ko.IssFF;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [22.5, -0.0862, 13.17, 1.0];
num_vars = 4;
lower_bd = [0,-2,0,0];
upper_bd = [50,1,30,2];
fit_fn = @(X) Iss_fitness(X,Iss_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 6);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = Iss(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value (steady-state current)
    fitted_Iss = A(:,66);
    
    % save the results
    file_path = sprintf('./results/GA_Iss_ko_%i.mat', i);
    rst(i,:) = [X, fval, fitted_Iss(end)];
    save(file_path, 'rst');
    disp(rst)
end


%% Iss - WT
% import data
K_wt = readtable("./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx");
Iss_wt = K_wt.Iss;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [22.5, -0.0862, 13.17, 1.0];
num_vars = 4;
lower_bd = [0,-2,0,0];
upper_bd = [50,1,30,2];
fit_fn = @(X) Iss_fitness(X,Iss_wt,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 6);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = Iss(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value (steady-state current)
    fitted_Iss = A(:,66);
    
    % save the results
    file_path = sprintf('./results/GA_Iss_wt_%i.mat', i);
    rst(i,:) = [X, fval, fitted_Iss(end)];
    save(file_path, 'rst');
    disp(rst)
end
