clc
close all
clear variables


%% Main
% parameters for voltage clamp
holding_p = -140; %mV
holding_t = 10; %msec
P1s = -130:10:50; %mV
P1_t = 80; % msec
P2 = -20; % mV
P2_t = 100; % msec

% run the model
t_container = cell(1,length(P1s));
S_container = cell(1,length(P1s));
A_container = cell(1,length(P1s));
for i=1:length(P1s)
    [t,S,A,C] = Rasmusson(holding_p,holding_t,P1s(i),P1_t,P2,P2_t);
    t_container{i} = t;
    S_container{i} = S;
    A_container{i} = A;
end

% plots
figure(1);
for i=1:length(P1s)
    hold on
    subplot(2,2,1);
    plot(t_container{i},A_container{i}(:,72));
    axis tight;
    xlabel('Time (msec)');
    ylabel('V (mV)');
    title('Two-Pulse Protocol');
    hold off
    
    % A40; Fast Na+ current
    hold on
    subplot(2,2,2);
    plot(t_container{i},A_container{i}(:,58));
    axis tight
    xlabel('time (msec)');
    ylabel('(pA/pF)');
    title('Fast Na+ Current I_{Na}');
    hold off
    
    % A3; Cass
    hold on
    subplot(2,2,3);
    plot(t_container{i},S_container{i}(:,3));
    axis tight
    xlabel('Time (msec)');
    ylabel('uM');
    title('[Ca^{2+}]_{ss}');
    hold off
    
    % A67; I_Kto,f
    hold on
    subplot(2,2,4);
    plot(t_container{i},A_container{i}(:,61));
    axis tight
    xlabel('Time (msec)');
    ylabel('pA/pF');
    title('I_{Kto,f}');
    hold off
end

% [t,S,A,C] = Ito(holding_p,holding_t,P1,P1_t,P2,P2_t,[30,30,13.5,33.5,33.5,33.5]);
% [t,S,A,C] = Rasmusson(holding_p,holding_t,P1,P1_t,P2,P2_t);

% visualization
% subplot(2,1,1)
% plot(t, A(:,66))
% subplot(2,1,2)
% plot(t, A(:,61))


%% transform .mat arrays to Excel files
file_names = dir('./results/*.mat');
for i=1:length(file_names)
    file_name = file_names(i);
    read_path = sprintf('./results/%s', file_name.name);
    load(read_path)
    
    [~, write_name, ~] = fileparts(file_name.name);
    write_path = sprintf('./results/%s.xlsx', write_name);
    xlswrite(write_path, rst)
end
