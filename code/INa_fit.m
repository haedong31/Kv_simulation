clc
close all
clear variables


%% INa
% import data
% na_gv_ko = readtable("./MGAT1_Data_tidy/JMCC/Nav Currents/Ina GV MGAT1KO Final.xlsx", ...
%     'Range', 'A1:T24', ...
%     'ReadVariableNames', true);
% na_gv_wt = readtable("./MGAT1_Data_tidy/JMCC/Nav Currents/Ina GV MGAT1KO Final.xlsx", ...
%     'Range', 'X1:AQ17', ...
%     'ReadVariableNames', true);
% na_ssi_ko = readtable("./MGAT1_Data_tidy/JMCC/Nav Currents/Ina SSI MGAT1KO Final.xlsx", ...
%     'Range', 'B1:U26', ...
%     'ReadVariableNames', true);
% na_ssi_wt = readtable("./MGAT1_Data_tidy/JMCC/Nav Currents/Ina SSI MGAT1KO Final.xlsx", ...
%     'Range', 'Z1:AS17', ...
%     'ReadVariableNames', true);

na_gv_ko = readtable("./Ina GV MGAT1KO Final.xlsx", ...
    'Range', 'A1:T24', ...
    'ReadVariableNames', true);
na_gv_wt = readtable("./Ina GV MGAT1KO Final.xlsx", ...
    'Range', 'X1:AQ17', ...
    'ReadVariableNames', true);
na_ssi_ko = readtable("./Ina SSI MGAT1KO Final.xlsx", ...
    'Range', 'B1:U26', ...
    'ReadVariableNames', true);
na_ssi_wt = readtable("./Ina SSI MGAT1KO Final.xlsx", ...
    'Range', 'Z1:AS17', ...
    'ReadVariableNames', true);

%% KO
SSA = table2array(na_gv_ko(:,7:20));
SSI = table2array(na_ssi_ko(:,5:20));
Gmax = 0.9;

num_vars = 25;
lower_bd = [-2.5, -2.5, -2.5, -7.5, -2.5, ...
            -12.5, 0, -2.5, -2.5, -2.5, ...
            2.0, -0.001, 2.0, -0.001, 2.0, ...
            -0.001, 2.0, 0.00001, -1.0, 0.00001, ...
            0.001, 3.2, 1.4, 1.4, 1.4];
upper_bd = [7.5, 7.5, 7.5, 2.5, 7.5, ...
            -2.5, 5.0, 7.5, 7.5, 7.5, ...
            12.0, 0.001, 12.0, 0.001, 12.0, ...
            0.001, 12.0, 0.009, 2.0, 0.009, ...
            0.9, 30, 14, 14, 14];
fit_fn = @(X) INa_fitness(X, SSA, SSI, Gmax);

% run GA
rst = zeros(10, 26);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % save the results
    file_path = sprintf('./GA_INa_ko_%i.mat', i);
    rst(i,:) = [X, fval];
    save(file_path, 'rst');
    disp(rst)
end


%% WT
SSA = table2array(na_gv_wt(:,7:20));
SSI = table2array(na_ssi_wt(:,5:20));
Gmax = 1.0;

num_vars = 25;
lower_bd = [-2.5, -2.5, -2.5, -7.5, -2.5, ...
            -12.5, 0, -2.5, -2.5, -2.5, ...
            2.0, -0.001, 2.0, -0.001, 2.0, ...
            -0.001, 2.0, 0.00001, -1.0, 0.00001, ...
            0.001, 3.2, 1.4, 1.4, 1.4];
upper_bd = [7.5, 7.5, 7.5, 2.5, 7.5, ...
            -2.5, 5.0, 7.5, 7.5, 7.5, ...
            12.0, 0.001, 12.0, 0.001, 12.0, ...
            0.001, 12.0, 0.009, 2.0, 0.009, ...
            0.9, 30, 14, 14, 14];
fit_fn = @(X) Iss_fitness(X, SSA, SSI, Gmax);

% run GA
rst = zeros(10, 26);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % save the results
    file_path = sprintf('./GA_INa_wt_%i.mat', i);
    rst(i,:) = [X, fval];
    save(file_path, 'rst');
    disp(rst)
end


% %% KO - SSA
% % clc
% % close all
% % clear variables
% 
% % SSA protocol
% holding_p = -100; % mV
% holding_t = 10; % ms
% P1 = -85:5:-20; % mV
% P1_t = 120 + holding_t; % ms
% P2 = -100; % mV
% P2_t = P1_t; % ms
% 
% % GA parameters
% X = [2.5, 2.5, 2.5, -2.5, 2.5, ...
%     -7.5, 2.5, 2.5, 2.5, 2.5, ...
%     7.0, 0.0084, 7.0, 0.0084, 7.0, ...
%     0.0084, 7.0, 1/1000.0, 1.0, 1/95000.0, ...
%     1/50.0, 16.6, 7.7, 7.7, 7.7];
% 
% % SSA parameters
% Gmax = 0.9;
% 
% % run simulation
% ts = cell(1, length(P1));
% As = cell(1, length(P1));
% for i=1:length(P1)
%     [t, ~, A, ~] = INa(holding_p,holding_t,P1(i),P1_t,P2,P2_t,X);
%     ts{i} = t;
%     As{i} = A;
% end
% 
% % extract peaks and Ernest potentials at peak
% peaks = zeros(1, length(P1));
% ENas = zeros(1, length(P1));
% for i=1:length(P1)
%     [peak, peak_idx] = min(As{i}(:,58));
%     peaks(i) = peak;
%     ENas(i) = As{i}(peak_idx,57);
% end
% 
% % G = I_peak/(V-ENa)
% % SSA = G/Gmax
% SSA = peaks./((P1 - ENas)*13);
% 
% 
% %% KO - SSI
% % clc
% % close all
% % clear variables
% 
% % SSI protocol
% holding_p = -100; % mV
% holding_t = 10; % ms
% P1 = -140:5:-65; % mV
% P1_t = 500 + holding_t; % ms
% P2 = -20; % mV
% P2_t = P1_t + 10; % ms
% 
% % GA parameters
% X = [2.5, 2.5, 2.5, -2.5, 2.5, ...
%     -7.5, 2.5, 2.5, 2.5, 2.5, ...
%     7.0, 0.0084, 7.0, 0.0084, 7.0, ...
%     0.0084, 7.0, 1/1000.0, 1.0, 1/95000.0, ...
%     1/50.0, 16.6, 7.7, 7.7, 7.7];
% 
% % run simulation
% ts = cell(1, length(P1));
% As = cell(1, length(P1));
% for i=1:length(P1)
%     [t, ~, A, ~] = INa(holding_p,holding_t,P1(i),P1_t,P2,P2_t,X);
%     ts{i} = t;
%     As{i} = A;
% end
% 
% % extract peaks
% peaks = zeros(1, length(P1));
% for i=1:length(P1)
%     [peak, peak_idx] = min(As{i}(:,58));
%     peaks(i) = peak;
% end
% 
% for i=1:length(P1)
%     hold on
%     plot(ts{i}, As{i}(:,58))
% end
% Imax = min(peaks);
% SSI = peaks/Imax;
