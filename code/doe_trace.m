clc
close all
clear variables


%% import relevant data
% DoE matrix
doe_mx = readtable('./FrF2_dgn_matrix.csv');

% trace
load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};


% real experiment data
K_data = readtable('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx');

Iss = K_data.IssFF;
Iss = nanmean(Iss);

cap = K_data.CapFF;
cap = nanmean(cap);

% normalize the trace
ds_Ktrace.I = ds_Ktrace.I ./ cap;


%% DoE
X0 = [-22.5, 7.7, -45.2, 5.7, 0.18064, 0.03577, 30.0, 0.3956, -0.06237, 30.0, ...
    0.000152, 13.5, 7.0, 0.067083, 33.5, 7.0, 0.00095, 33.5, 7.0, 0.051335, ...
    22.5, 7.7, 26.5, 7.7, 45.2, 5.7, 0.493, -0.0629, 2.058, -170.0, 45.2, 5.7, ...
    22.5, 7.7, 26.5, 7.7, 45.2, 5.7, 0.493, -0.0629, 2.058, -170.0, 45.2, 5.7];

IKsum = Kv_anal(t, 50, X0);
plot(t, IKsum)
num_vars = length(X0);
low = zeros(1, num_vars);
high = zeros(1, num_vars);
for i=1:length(X0)
    low(i) = X0(i) - 3*abs(X0(i));
    high(i) = X0(i) + 3*abs(X0(i));
end

num_runs = size(doe_mx);
num_runs = num_runs(1);
res = zeros(num_runs, 1); 
for i=1:num_runs
    treat = table2array(doe_mx(i, :));
    param = zeros(1, num_vars);
    for j=1:num_vars
        if treat(j) == 1
            param(j) = high(j);
        else
            param(j) = low(j);
        end
    end
    
    t = ds_Ktrace.time;
    IKsum = Kv_anal(t, 50, param);
    
    [peak, peak_idx] = max(abs(IKsum));
    if peak ==  Inf
        IKsum(peak_idx:end) = IKsum(peak_idx-1);
    else
        IKsum(peak_idx:end) = IKsum(peak_idx:end) + Iss;
    end
    
    trace_sim = dtw(IKsum, ds_Ktrace.I);
    res(i) = trace_sim;
end
