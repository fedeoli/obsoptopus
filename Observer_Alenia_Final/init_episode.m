%% init episode for RL
%%% workspace %%%
% clear
if exist('RL','var')
    keep RL s current_reward test_flag init_flag
    setup.load_mem = 0;
else
    clear
end
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DYNOPT INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% constant stuff %%%%%

% magnetometers displacement 
setup.RPYbetweenMagnetometers = 1*[0,0,90]*pi/180;

% state duration and main setup
setup.w = 5;
setup.Nts = 3;
setup.d = 0;
setup.u_amp = 0;
setup.lowpass_pwm = 0;
% short period setup
setup.T_duration = 5;
setup.t_start = 0;
setup.Tend = setup.t_start + setup.T_duration;
setup.Ts = 1e0;

% control
setup.control = 1;
setup.T = 100;
setup.u_freq = 1/setup.T;

% plot and graphics
setup.plot = 0;
setup.print = 0;

%%% integration %%%
setup.integration_pos = 1;
setup.integration_att = 1;

%%%%% OBSERVER %%%%%
setup.generate_plant = 1;
setup.generate_orbit = 1;
setup.ObserverOn = 1;
setup.Observer = 'OPT';
setup.OptimisationOn = 1;
setup.simulationModel = 1;
setup.identify = 1;
setup.check = 0;
setup.fault_sim = 0;
setup.flush_buffer = 0;
setup.always_opt = 1;
setup.optimise_input = 0;
setup.input_tuning = 0;

% conditions to consider the optimisation result
setup.Jdot_thresh = 9e-1;
setup.blue_flag = 0;
% built in/gradient optimisation conditions
setup.J_thresh = [1e-10, 1e3];
setup.max_iter = 20;
setup.maxFcount = Inf;
setup.safety_density = 1;

%%%% HYSTERESIS %%%%
% the optimisation is run if COND > THRESH (set to 0 for a nocare condition)
setup.adaptive = 1;
setup.dJ_2 = setup.adaptive*1e-2;
setup.dJ_1 = setup.adaptive*5e-3;

% optimisation
setup.fcon_flag = 0;
if setup.fcon_flag == 1
    setup.fmin = @fmincon;
else
    setup.fmin = @fminsearch;
end

%%%%% INTEGRATION %%%%%
setup.forward = 1;
%%%% initialise DynOpt %%%%     
setup.load_mem = 0;
%%%% setup as it was the first algorithm iteration (MDP state init) %%%
setup.RL = 1;

%%%%% MODEL %%%%%
setup.model = 'satellite';

%% episode
if init_flag == 1
    % iteration
    RL.S.i = 1;

    % action
    RL.S.A(:,1) = RL.A.domain_A(:,1);

    % initial attitude - true
    RL.S.satellites_attitude_true = (RL.E.domain_status(:,2)-RL.E.domain_status(:,1)).*rand(RL.E.dimState,1) + RL.E.domain_status(:,1);

    % initial attitude - est
    setup.noise_enable = 1;

    % error init definition
    %%%%% NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    setup.ErrorOnEA = 1;
    setup.rand_init = 0;
    % params estimation
    setup.nbias = 0;
    setup.bias_dyn = 0;
    setup.bias_enable = 0;
    setup.bias_mag_enable = 0;
    setup.optimise_params = 1;
    setup.nparams = 3;
    setup.inertia = 1;
    % bias
    if setup.bias_enable 
        setup.bias = 1*ones(3,1)*5*pi/180; 
    else
        setup.bias = zeros(3,1);
    end
    if setup.bias_mag_enable 
        setup.bias_mag = 1*abs(1e-4*randn(6,1)+1e-3);
    else
        setup.bias_mag = zeros(6,1);
    end
    setup.bias_tot = setup.bias;

    % measures
    setup.EulerAngleNoiseOnMag = 0*setup.noise_enable*1e-2;
    setup.noise_amp = 1*setup.noise_enable*[1*1e-4*ones(3,1); 1*1e-4*ones(6,1)];
    %%%% STATE %%%%
    if setup.integration_pos == 1 && setup.integration_att == 1
        % POS + ATT
        setup.init_error_amp = 1*setup.noise_enable*[0*ones(3,1); 0*ones(3,1); 20*ones(4,1); 5*ones(3,1)];
    elseif setup.integration_pos == 1 && setup.integration_att == 0
        % POS
        setup.init_error_amp = 1*setup.noise_enable*[1e1*ones(3,1); 5e-1*ones(3,1)];
    else
        % ATT (errors in percentage)
        setup.init_error_amp = 1*setup.noise_enable*[40*ones(4,1); 10*ones(3,1)];
    end
    %%%% PARAMS %%%%
    setup.ParamsAllPos = 1;
    setup.init_param_error_amp = 1*setup.noise_enable*ones(1,setup.nparams)*50;

    % init state
    if setup.noise_enable
        RL.S.S0 = (RL.E.domain_status(:,2)-RL.E.domain_status(:,1)).*rand(RL.E.dimState,1) + RL.E.domain_status(:,1);
    else
        RL.S.S0 = RL.S.satellites_attitude_true;
    end 

    % target attitude
    RL.S.T0 = (RL.E.domain_target(:,2)-RL.E.domain_target(:,1)).*rand(RL.E.dimTarget,1) + RL.E.domain_target(:,1);
    setup.RL_data = RL;
    setup.target_attitude = RL.S.T0;
    
    %%%%%%%%%%%% random environment definition %%%%%%%%%%%%%%%%%%%%%%%%
    %%% ORBIT %%%
    ecc = (RL.E.domain_ecc(2)-RL.E.domain_ecc(1)).*rand(1) + RL.E.domain_ecc(1);
    inclination = (RL.E.domain_i(2)-RL.E.domain_i(1)).*rand(1) + RL.E.domain_i(1);
    om = (RL.E.domain_om(2)-RL.E.domain_om(1)).*rand(1) + RL.E.domain_om(1);
    RAAN = (RL.E.domain_RAAN(2)-RL.E.domain_RAAN(1)).*rand(1) + RL.E.domain_RAAN(1);
    f0 = (RL.E.domain_f0(2)-RL.E.domain_f0(1)).*rand(1) + RL.E.domain_f0(1);
    T = (RL.E.domain_T(2)-RL.E.domain_T(1)).*rand(1) + RL.E.domain_T(1);
    RL.S.orbit(:,RL.S.i) = [ecc; inclination; om; RAAN; f0; T];
    setup.orbit = RL.S.orbit(:,RL.S.i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% get first nu and iner_ECI %%%%
    setup.Tend = 2;
    [DynOpt, ~] = ObsOpt_RL_v1_fun(setup);   
else
    % orbit
    RL.S.orbit(:,RL.S.i) = RL.S.orbit(:,RL.S.i-1);
    setup.orbit = RL.S.orbit(:,RL.S.i);
end