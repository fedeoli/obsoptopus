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

% state
DynOpt.integration_pos = 1;
DynOpt.integration_att = 1;
params.Control = 1;

% observer
params.observed_state = [5 6 7];

if DynOpt.integration_pos == 1 & DynOpt.integration_att == 0
    DynOpt.StateDim = 6*params.Nagents;
    DynOpt.StateDim_single_agent = 6;
    DynOpt.init_state = satellites_iner_ECI;
    DynOpt.dim_out = length(params.observed_state);
    DynOpt.scale_factor = 1e1*ones(1,DynOpt.dim_out);
elseif DynOpt.integration_pos == 0 & DynOpt.integration_att == 1
    DynOpt.StateDim = 7*params.Nagents;
    DynOpt.StateDim_single_agent = 7;
    DynOpt.init_state = satellites_attitude;
    DynOpt.dim_out = length(params.observed_state)+3*2;
    DynOpt.scale_factor = ones(1,DynOpt.dim_out);
elseif DynOpt.integration_pos == 1 & DynOpt.integration_att == 1
    DynOpt.StateDim = 13*params.Nagents;
    DynOpt.StateDim_single_agent = 13;
    DynOpt.init_state = [satellites_iner_ECI; satellites_attitude];
    DynOpt.dim_out = length(params.observed_state)+3*2;
    DynOpt.scale_factor = ones(1,DynOpt.dim_out);
end

DynOpt.param_estimate = [params.sat(1).I(1,1)];
params.param_estimate = DynOpt.param_estimate;


% integration time
DynOpt.Niter = length(time);
DynOpt.time = time;
DynOpt.Ts = 1;
DynOpt.Tstart = time(1);
DynOpt.Tend = time(end);

DynOpt.magnetoBias = zeros(3,1);
DynOpt.gyroBias = zeros(3,1);
DynOpt.RPYbetweenMagSensors = [0; 0; pi/2];

% params
DynOpt.params = params;
DynOpt.params.input_flag = params.Control;
params.eps_coef = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of output arrays
OutputInitialization_V2_3;

if params.Observer    
    [ObserverTest, Agent] = SetObserver_v1_1(1, length(DynOpt.time), satellites_iner_ECI, satellites_attitude, params); 
end

% check attitude motion 
params.Attitude = DynOpt.integration_att;

% check noise on attitude
if DynOpt.noise_enable == 0
   ObserverTest.AttitudeZeroErrors = 1; 
else
   ObserverTest.AttitudeZeroErrors = 0; 
end



