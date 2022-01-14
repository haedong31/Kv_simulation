% pilot parameter estimation of the Boltzmann equation
% steady-state activation and inactivation of Na channel

clear variables
close all
clc


%% import the data
ssa=readtable('./Nav_WT_ggmax.csv');
voltages=ssa.voltage;
ss_acts=ssa.mean_ggmax;

% define a Boltzmann function
boltzmann = @(v, theta1, theta2) 1/(1+exp(-(v-theta1)/theta2));

% algorihm parameters
iters=100000;
num_data=size(ssa,1);
psig=1;
theta1=zeros(1,iters);
theta2=zeros(1,iters);
errors=zeros(1,iters);

% initialization
theta1(1)=0;
theta2(1)=0;

for i=1:iters
    % a new proposed parameters
    new_theta1=random('normal',-40,psig);
    new_theta2=random('normal',10,psig);
    
    % evaluation
    old_errs=zeros(1,num_data);
    new_errs=zeros(1,num_data);
    for j=1:num_data
        old_errs(j)=(boltzmann(voltages(j),theta1(i),theta2(i))-ss_acts(j))^2;
        new_errs(j)=(boltzmann(voltages(j),new_theta1,new_theta2)-ss_acts(j))^2;
    end
    errors(i)=sum(new_errs);
    
    % ratio of normal likelihood function
    ratio=exp(-sum(new_errs)+sum(old_errs));
    
    if rand<ratio
        theta1(i+1)=new_theta1;
        theta2(i+1)=new_theta2;
    else
        theta1(i+1)=theta1(i);
        theta2(i+1)=theta2(i);
    end
end

% plot with the simulated paramters
[min_err, min_idx]=min(errors);
vv=min(voltages):0.01:max(voltages);
fit_vals=zeros(1,length(vv));
for k=1:length(vv)
    fit_vals(k)=boltzmann(vv(k),theta1(min_idx),theta2(min_idx));
    % fit_vals(k)=boltzmann(vv(k),22.5,7);
end    
plot(vv,fit_vals,'Color','b','LineWidth',2)

% plot the real data
hold on
plot(voltages,ss_acts,'-ro','Color','r','LineWidth',2)
hold off
legend('Simulated', 'Real')
xlabel('Voltage (mV)')
ylabel('SS Activation (WT)')
