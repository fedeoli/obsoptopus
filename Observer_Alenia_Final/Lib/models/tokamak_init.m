%% PLANT model and data

global DynOpt params

Ki = 1.5; %good 1.5*2;
Kp = 0.22; %good 0.22*30;
Kd = 0.0; %good 1.0;
beta = 140;%good 30;
u0 = 1.0;

% model simulation
params.c0 = 45;
params.c1 = 20;
params.gamma0 = 10;
params.gamma1 = 0.01;
params.beta = beta;
params.StateDim = 4;
params.Zh = 1.5;
params.Zf = 1.3;
params.beta = beta;
params.u0 = u0;
params.gainPHD = 1/(1.5e-4);

params.StateDim = 4;
params.observed_state = 3;

params.input_flag = 1;

if params.input_flag == 0
    params.Ki = 0;
    params.Kp = 0;
    params.Kd = 0;
    params.u0 = 0;
    params.gainPHD = 0;
else
    params.Ki = Ki*1;
    params.Kp = Kp*1;
    params.Kd = Kd*1;
    params.u0 = u0;
    params.gainPHD = 1/(1.5e-4);
    
    % useless in this model - just to avoid errors    
    params.U = zeros(1,DynOpt.Niter);
end

DynOpt.params = params;