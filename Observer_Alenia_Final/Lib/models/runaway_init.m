%% PLANT model and data
% model simulation
% plant data

global DynOpt params

% electric charge
params.Q = 1;

% spontaneous emission
params.S = 0.01;

% parameters
params.gamma = 0.5;
params.gamma1 = 4;
params.ni = 0.5;
params.Wt = 0.1;

% ringing
params.wq = 1;
params.chi = 1;

% eps_coef
params.eps_coef = 1;

% initial condition
% params.T0 = 9.255;
% params.W0 = 0.02176;
params.T0 = 0;
params.W0 = 0.001;

params.StateDim = 2;

if DynOpt.simulationModel == 1
    params.observed_state = 2;
else
    params.observed_state = 2;
end

params.input_flag = 0;

if params.input_flag == 0
    params.U = zeros(1,DynOpt.Niter);
else
    params.U = zeros(1,DynOpt.Niter);
end

DynOpt.params = params;