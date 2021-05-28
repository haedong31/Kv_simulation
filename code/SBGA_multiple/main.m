clc
close all
clear variables

% read data
wt_data = table2array(readtable('4.5s-avg-wt.csv'));
[~, num_volts] = size(wt_data);
num_volts = num_volts - 1;

% calibration arguments
tol = [0.1, 1.0];

N0 = 30;
N1 = 6;
N2 = 4;

% run calibration
par_to_wt = ito_calibration(Ato_wt, tto_wt, tol, N0, N1, N2);
par_kslow_wt = ikslow_calibration(Akslow_wt, tkslow_wt, tol, N0, N1, N2);

par_to_ko = ito_calibration(Ato_ko, tto_ko, tol, N0, N1, N2);
par_kslow_ko = ikslow_calibration(Akslow_ko, tkslow_ko, tol, N0, N1, N2);
