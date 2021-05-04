%% PLANT model and data
% model simulation
% plant data

global DynOpt params

params.Lt = 1;
params.M = 1;
params.g = -9.81;
params.Fw = 0.5;
params.Fc = 0.5;
params.If = 5;
params.Iw = 3;
params.Km = 1;

params.StateDim = 4;
params.observed_state = 4;

params.input_flag = 1;
params.Tm = params.Km*sin(20*DynOpt.time);% + params.Km*sin(15*DynOpt.time);
params.U = params.Tm;

DynOpt.params = params;