clc
close all
clear variables


%% code parameters
% voltage clamp protocol parameters
holding_p = -70;
holding_t = 4.5; 
P1 = 50;
P1_t = 25*1000;
P2 = 50; 
P2_t = P1_t;

init_size = 20;


%% initialization
% original values in the Rasmusson: IKslow
X_slow = [22.5000, 2.05800, 45.2000, 1200.00, 45.2000, 0.493000, 170.000];
X_to = [30, 30, 13.5, 33.5, 33.5, 33.5];

num_params = length(X_slow) + length(X_to);
lower_bd = zeros(1, num_params);
upper_bd = zeros(1, num_params);
for i=1:(length(X_slow) + length(X_to))
    lower_bd(i) = X(i) - X(i)*2;
    upper_bd(i) = X(i) + X(i)*2;
end

% sampling
init_params = zeros(init_size, num_params);
for i=1:num_params
    unif = makedist('Uniform', 'lower',lower_bd(i), 'upper',upper_bd(i));
    init_params(:,i) = random(unif, init_size, 1);
end

y = ds_Ktrace_ko;
y_len = length(ds_Ktrace_ko.KO);
hat = table(t, A(:,65));
[peak, peak_idx] = max(hat.Var2);
sub_hat = hat(peak_idx:end,:);
ds_sub_hat = downsample(sub_hat, floor(length(t)/y_len));
ds_sub_hat = ds_sub_hat(1:y_len,:);
sse = sum((y.KO - ds_sub_hat.Var2).^2);
