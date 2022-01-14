clear variables
close all
clc


%% import the AP Data WT
opts = delimitedTextImportOptions("NumVariables", 27);

% Specify range and delimiter
opts.DataLines = [2, 20];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VAR1", "VAR2", "MemPot", "VarName4", "Hold", "Peak", "TTP", "MaxSlope", "TMS", "AP25", "AP50", "AP75", "AP90", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27"];
opts.SelectedVariableNames = ["VAR1", "VAR2", "MemPot", "VarName4", "Hold", "Peak", "TTP", "MaxSlope", "TMS", "AP25", "AP50", "AP75", "AP90"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [2, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
AP_WT = readtable("./MGAT1_Data_tidy/JMCC/APs/AP Parameters FF.CSV", opts);
clear opts


%% import the Ca2+ Data
opts = delimitedTextImportOptions("NumVariables", 23);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Date", "CellCONTROL", "Widthms", "HeightV", "TimeToPeakms", "MaxSlopeVs", "TMaxSlopems", "Width10ms", "Width50ms", "Width70ms", "Taus", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23"];
opts.SelectedVariableNames = ["Date", "CellCONTROL", "Widthms", "HeightV", "TimeToPeakms", "MaxSlopeVs", "TMaxSlopems", "Width10ms", "Width50ms", "Width70ms", "Taus"];
opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Ca_WT = readtable("./MGAT1_Data_tidy/JMCC/Ca Imaging 37 Degrees/Ca Imaging MGAT1KO Final.CSV", opts);
clear opts
