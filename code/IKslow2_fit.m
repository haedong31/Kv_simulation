clc
close all
clear variables


%% KO
K_ko = readtable("./potassium-KO.xlsx");
IKslow2_ko = K_ko.A2FF;

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
opts = optimoptions('ga', 'InitialPopulationRange', [-0.407070366; 7.284902402; 6.417898622; 3.54893446; 2.095174939; 1.542657627; 0.483600948], ...
    'PlotFcn', {@gaplotbestf,@gaplotdistance});
fit_fn = @(X) IKslow2_fitness(X,IKslow2_ko,holding_p,holding_t,P1,P1_t,P2,P2_t,opts);

% run GA
rst = zeros(10, 9);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[]);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow2(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow2 = A(:,65);
    
    % save the results
    file_path = sprintf('./GA_IKslow2_ko_%i.mat', i);
    rst(i,:) = [X, fval, max(fitted_IKslow2)];
    save(file_path, 'rst');
    disp(rst)
end


%% WT
K_wt = readtable("./potassium-WT.xlsx");
IKslow2_wt = K_wt.A2;

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
opts = optimoptions('ga', 'InitialPopulationRange', [-7.794253319, 1.250916659; 5.486307644; 13.47305609; 6.419422564; -2.034281757; -2.974514692], ...
    'PlotFcn', {@gaplotbestf,@gaplotdistance});
fit_fn = @(X) IKslow2_fitness(X,IKslow2_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 9);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],opts);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow2(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow2 = A(:,65);
    
    % save the results
    file_path = sprintf('./GA_IKslow2_wt_%i.mat', i);
    rst(i,:) = [X, fval, max(fitted_IKslow2)];
    save(file_path, 'rst');
    disp(rst)
end
