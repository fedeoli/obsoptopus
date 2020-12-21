%% PLANT model and data
% model simulation
% plant data

global DynOpt params 
if DynOpt.check == 0
    % electric charge
    params.Q = 1;

    % spontaneous emission
    params.S = 0;

    % parameters
    if struct.identify == 1
        params.gamma = 0.5;
        params.gamma1 = 2.5;
        params.ni = 0.5;
        params.Wt = 0.1;
    else
        params.gamma = struct.gamma;
        params.gamma1 = struct.gamma1;
        params.ni = 0.5;
        params.Wt = struct.Wt;
    end

    % eps_coef
    params.eps_coef = 1;

    % initial condition
    params.T0 = 0;
    params.W0 = 0.001;

    %%%%%%%%% WORKSPACE %%%%%%%%%%
    % params.T0 = 9.255;
    % params.W0 = 0.02176;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
else
    clear params
    save current
    load simulation/runaway/data_manual/default_510_con
    params = DynOpt.params_estimate_final;
    params.T0 = DynOpt.OptXstory(1,1);
    params.W0 = DynOpt.OptXstory(2,1);
    params.eps_coef = 1;
    keep params
    load current
end

DynOpt.params = params;