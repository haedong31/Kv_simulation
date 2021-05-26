clc
close all
clear variables

wt_data = readtable('Ik-4.5s-WT.csv');
ko_data = readtable('Ik-4.5s-KO.csv');

Ato_wt = wt_data.peak_Ito;
tto_wt = wt_data.tau_Ito;
Akslow_wt = wt_data.peak_IKslow;
tkslow_wt = wt_data.tau_IKslow;

Ato_ko = ko_data.peak_Ito;
tto_ko = ko_data.tau_Ito;
Akslow_ko = ko_data.peak_IKslow;
tkslow_ko = ko_data.tau_IKslow;

tol = [0.1, 1.0];

N0 = 30;
N1 = 6;
N2 = 4;

par_to_wt = ito_calibration(Ato_wt, tto_wt, tol, N0, N1, N2);
par_kslow_wt = ikslow_calibration(Akslow_wt, tkslow_wt, tol, N0, N1, N2);

par_to_ko = ito_calibration(Ato_ko, tto_ko, tol, N0, N1, N2);
par_kslow_ko = ikslow_calibration(Akslow_ko, tkslow_ko, tol, N0, N1, N2);
