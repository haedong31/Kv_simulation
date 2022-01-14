clc
close all
clear variables

% initialization
voltages=-120:1:-20;
ggmax_control = zeros(1,length(voltages));
ggmax_mgat1 = zeros(1,length(voltages));
iimax_control = zeros(1,length(voltages));
iimax_mgat1 = zeros(1,length(voltages));

% inactivation parameters
Vi_control = -82.4;
Ki_control = -5.5;
Vi_mgat1 = -78.7;
Ki_mgat1 = -6.7;

% activation parameters
Va_control = -47.6;
Ka_control = 5.2;
Va_mgat1 = -42.9;
Ka_mgat1 = 5.8;

for n=1:length(voltages)
    iimax_control(n)=1/(1+exp(-(voltages(n)-Vi_control)/Ki_control));
    iimax_mgat1(n)=1/(1+exp(-(voltages(n)-Vi_mgat1)/Ki_mgat1));
    ggmax_control(n)=1/(1+exp(-(voltages(n)-Va_control)/Ka_control));
    ggmax_mgat1(n)=1/(1+exp(-(voltages(n)-Va_mgat1)/Ka_mgat1));
end

plot(voltages, ggmax_control, 'k','LineWidth', 2)
hold on 
plot(voltages, ggmax_mgat1,'r','LineWidth',2)
plot(voltages, iimax_control, 'k', 'LineWidth', 2)
plot(voltages, iimax_mgat1,'r','LineWidth', 2)
hold off
legend('Control Group', 'MGAT1KO')
xlabel('Voltage (mV)')
ylabel('I/I_{max}, G/G_{max}')
