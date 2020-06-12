clc
close all
clear variables

y = [24.8, 105.2];
[bamp, btau, best_chroms] = Ito_AGA(5, y, 30, 6, 4);
