%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     ANALYSIS SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DynOpt params ObserverTest Agent

Scenario_K_ORB_A
params.Observer = DynOpt.ObserverOn;
DynOpt.ObserverTest = ObserverTest;
DynOpt.Agent = Agent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               DYNOPT INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% state
DynOpt.integration_pos = 1;
DynOpt.integration_att = 1;
DynOpt.control = struct.control;

if DynOpt.control == 1
    params.Control = 1;
end

%% observer
params.observed_state = [5 6 7];
DynOpt.n_sensor = 0;

if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
    DynOpt.StateDim = 6*params.Nagents;
    DynOpt.StateDim_single_agent = 6;
    DynOpt.init_state = satellites_iner_ECI;
    DynOpt.dim_out = length(params.observed_state);
elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
    DynOpt.StateDim = 7*params.Nagents;
    DynOpt.StateDim_single_agent = 7;
    DynOpt.init_state = satellites_attitude;
    DynOpt.dim_out = length(params.observed_state)+3*DynOpt.n_sensor;
elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
    DynOpt.StateDim = 13*params.Nagents;
    DynOpt.StateDim_single_agent = 13;
    DynOpt.init_state = [satellites_iner_ECI; satellites_attitude];
    DynOpt.dim_out = length(params.observed_state)+3*2;
end

%% integration time
DynOpt.Niter = length(time);
DynOpt.time = time;
DynOpt.Ts = 1;
DynOpt.Tstart = time(1);
DynOpt.Tend = time(end);
DynOpt.tspan = [1, 1+DynOpt.Ts];

%% model
params.eps_coef = 1;
params.bias = 10*pi/180; %from [deg/sec] to [rad/sec]

%%%%%% GYROSCOPE NOISE MODEL %%%%
params.RW_var = 3e-4;
params.RW_mean = 0;
DynOpt.RW_mem = zeros(DynOpt.Niter,1);
params.mu_bias = 0;

%% PARAMETERS %%
% simulation fault
DynOpt.fault_sim = struct.fault_sim;

% uncomment to estimate inertia
% DynOpt.param_estimate = [params.sat(1).I(1,1), params.sat(1).I(2,2), params.sat(1).I(3,3)];

% uncomment to estimate gyro bias
DynOpt.param_estimate = params.bias;
DynOpt.bias_dyn = 1;

% uncomment to estimate mass
% DynOpt.param_estimate = params.sat(1).M;
%%%%%%%%%%%%%%%%%

% old vars
% DynOpt.magnetoBias = zeros(3,1);
% DynOpt.gyroBias = zeros(3,1);
DynOpt.RPYbetweenMagSensors = [0; 0; pi/2];

%% update and init
params.param_estimate = DynOpt.param_estimate;
DynOpt.param_estimate_init = DynOpt.param_estimate;

% final params
DynOpt.params = params;

% control
DynOpt.params.input_flag = params.Control;

% nparams
DynOpt.nparams = length(DynOpt.param_estimate);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of output arrays
OutputInitialization_V2_3;

if params.Observer    
    [ObserverTest, Agent] = SetObserver_v1_1(1, length(DynOpt.time), satellites_iner_ECI, satellites_attitude, params); 
end

% check attitude motion 
params.Attitude = DynOpt.integration_att;
params.DesiredAttitude = zeros(3,1);

% check noise on attitude
if DynOpt.noise_enable == 0
   ObserverTest.AttitudeZeroErrors = 1; 
else
   ObserverTest.AttitudeZeroErrors = 0; 
end




