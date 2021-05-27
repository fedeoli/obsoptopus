function [DynOpt, params,satellites_iner_ECI,satellites_attitude] = satellite_init_function(DynOpt, params, struct)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     ANALYSIS SETTINGS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if DynOpt.generate_orbit == 0
        [params, satellites_iner_ECI, satellites_attitude] = Scenario_K_ORB_A_function(params,DynOpt,struct);
    else
        [params, satellites_iner_ECI, satellites_attitude] = Scenario_ObsOpt_function_v2(params,DynOpt,struct);
    end
    params.Observer = DynOpt.ObserverOn;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               DYNOPT INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% STATE AND OBSERVER %%%%
    DynOpt.integration_pos = struct.integration_pos;
    DynOpt.integration_att = struct.integration_att;
    DynOpt.control = struct.control;
    params.Control = DynOpt.control;

    %%% observation data %%%
    params.observed_state = DynOpt.integration_pos*6 + [5 6 7];
    DynOpt.n_sensor = 3;

    %%% state data %%%%
    if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
        DynOpt.StateDim = 6*params.Nagents;
        DynOpt.StateDim_single_agent = 6;
        DynOpt.init_state = satellites_iner_ECI;
    elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
        DynOpt.StateDim = 7*params.Nagents;
        DynOpt.StateDim_single_agent = 7;
        DynOpt.init_state = satellites_attitude;
    elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
        DynOpt.StateDim = 13*params.Nagents;
        DynOpt.StateDim_single_agent = 13;
        DynOpt.init_state = [satellites_iner_ECI; satellites_attitude];
    end

    %%% INTEGRATION SETUP %%%%
    DynOpt.Niter = length(params.time);
    DynOpt.time = params.time;
    DynOpt.Ts = struct.Ts;
    DynOpt.Tstart = params.time(1);
    DynOpt.Tend = params.time(end);
    DynOpt.tspan = [1, 1+DynOpt.Ts];

    %%% MODEL PARAMETERS %%%%
    params.eps_coef = 1;
    params.bias = struct.bias_tot;
    params.MagnetoBias = struct.bias_mag;
    DynOpt.RPYbetweenMagSensors = struct.RPYbetweenMagnetometers;
    DynOpt.EulerAngleNoiseOnMag = struct.EulerAngleNoiseOnMag;
    DynOpt.fault_sim = struct.fault_sim;

    %%%%%% GYROSCOPE NOISE MODEL %%%%
    % nparams
    DynOpt.nparams = struct.nparams;
    params.RW_var = 1*3e-3;
    params.RW_var_mag = 1*3e-7;
    params.RW_mean = 0;
    params.RW_mean_mag = 0;
    DynOpt.RW_mem = zeros(DynOpt.Niter,DynOpt.nparams);
    params.mu_bias = 0;

    %%% PARAMETERS TO RECONSTRUCT %%%
    % uncomment to estimate inertia
    if DynOpt.inertia 
        % 1 inertia
%         inertia = params.sat(1).I(1,1);
        
        % 3 inertias
        inertia = [params.sat(1).I(1,1); params.sat(1).I(2,2); params.sat(1).I(3,3)];
    else
        inertia = [];
    end

    % uncomment to estimate gyro bias
    if DynOpt.nbias > 0
        bias = params.bias;
    else
        bias = [];
    end


    % final estimation setup
    DynOpt.param_estimate = [bias; inertia];

    %% UPDATE AND INIT %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.param_estimate = DynOpt.param_estimate;
    DynOpt.param_estimate_init = DynOpt.param_estimate;

    % final params
    DynOpt.params = params;

    % control
    DynOpt.params.input_flag = params.Control;

    %%%%%%%%% Initialization scripts %%%%%%%%%%
    params = OutputInitialization_V2_3_function(params,satellites_iner_ECI,satellites_attitude);
    if params.Observer    
        [ObserverTest, Agent] = SetObserver_v1_1(1, length(DynOpt.time), satellites_iner_ECI, satellites_attitude, params); 
        DynOpt.Agent = Agent;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% set orbital position %%%%%
    params.SatellitesCoordinates = satellites_iner_ECI;
    params.SatellitesAttitude = satellites_attitude;

    %%%%%%%% init magnetic field %%%%%%%%
    myutc = [2019 12 15 10 20 36]; %CHANGE THIS...??
    LatLongAlt = eci2lla(params.SatellitesCoordinates(1:3)'*1E3,myutc); %converto from ECI to latitude, longitude,  altitude
    [mag_field_vector,~,~,~,~] = igrfmagm(max(1000,min(LatLongAlt(3),6E5)),LatLongAlt(1),LatLongAlt(2),decyear(2019,12,15),13); %mag_field_vector is in nanotesla, by IGRF11-12
    DynOpt.mag_field_vector = mag_field_vector;

    %%%%% TEST ORC %%%%%
    % ChiefAttitude_VTEST;
    % params.SatellitesAttitude = satellites_attitude;


    %%%%% check attitude motion %%%%%
    params.Attitude = DynOpt.integration_att;
    params.DesiredAttitude = zeros(3,1);

    %%%%% check noise on attitude %%%%
    if DynOpt.noise_enable == 0
       ObserverTest.AttitudeZeroErrors = 1; 
    else
       ObserverTest.AttitudeZeroErrors = 0; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save ObserverTest
    DynOpt.ObserverTest = ObserverTest;
    
    % dim out
    DynOpt.dim_out  = 9;
end




