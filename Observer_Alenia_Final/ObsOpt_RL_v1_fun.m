function [DynOpt, params] = ObsOpt_RL_v1_fun(struct)

%% Init Section
% close all
if struct.print 
    clc
end

% new random seed
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%% global variables %%%
DynOpt = [];
params = [];

% generate plant
DynOpt.generate_plant = struct.generate_plant;

% Inter-sampling Time, number of samples within intersampled data
DynOpt.Nts = struct.Nts; 

%number of inter-sampled data in a window: the total number of sampled data in a window are (DynOpt.w-1)*DynOpt.Nts+1
DynOpt.w = struct.w;

% observer optimisation options
DynOpt.fcon_flag = struct.fcon_flag;
DynOpt.max_iter = struct.max_iter;
DynOpt.identify = struct.identify;

%%% ENABLE INPUT TUNING MODE %%%
DynOpt.input_tuning = struct.input_tuning;

% identify
DynOpt.identify = struct.identify;
DynOpt.bias_dyn = struct.bias_dyn;
DynOpt.bias_enable = struct.bias_enable;
DynOpt.bias_mag_enable = struct.bias_mag_enable;

%%% init flags and setup %%%
% model name
DynOpt.modelname = struct.model;

% observer setup
DynOpt.ObserverOn = struct.ObserverOn;
DynOpt.OptimisationOn = struct.OptimisationOn;

% starts finding the state/parameters backward, -1, (actual time t_k) or Forward, 1 (i.e. at time  t_{k-(w-1)*Nts}
DynOpt.ForwardOptimization = struct.forward; 

% measure noise
DynOpt.noise_enable = struct.noise_enable;
DynOpt.measure_amp = struct.noise_amp;
DynOpt.EulerAngleNoiseOnMag = struct.EulerAngleNoiseOnMag;

% simulate the plant
DynOpt.simulationModel = struct.simulationModel;

% wrap function
DynOpt.wrap = @wrapToPi;

% print option
DynOpt.print = struct.print;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLANT model and synthetic data
%%% INPUT SETUP%%%
DynOpt.u_amp = struct.u_amp;
DynOpt.u_freq = struct.u_freq;
DynOpt.d = struct.d;
DynOpt.target_attitude = struct.target_attitude;
DynOpt.set_input = @set_input_v5;
DynOpt.switch_pwm = 1;
DynOpt.t_lowpass = 1;
DynOpt.lowpass_pwm = struct.lowpass_pwm;

%%%% OUTPUT SETUP %%%%
DynOpt.nbias = struct.nbias;
DynOpt.inertia = struct.inertia;

%%%%%%%%%%% INIT SECTION %%%%%%%%%% 
%%% Orbit generation data %%%
DynOpt.generate_orbit = struct.generate_orbit;
DynOpt.domain_ecc = [1e-4; 2e-4];
DynOpt.domain_i = [0;pi/2];
DynOpt.domain_om = [0;pi/2];
DynOpt.domain_RAAN = [0;2*pi];
DynOpt.domain_f0 = [0;2*pi];
DynOpt.domain_T = [5e3;7e3];
ecc = (DynOpt.domain_ecc(2)-DynOpt.domain_ecc(1)).*rand(1) + DynOpt.domain_ecc(1);
inclination = (DynOpt.domain_i(2)-DynOpt.domain_i(1)).*rand(1) + DynOpt.domain_i(1);
om = (DynOpt.domain_om(2)-DynOpt.domain_om(1)).*rand(1) + DynOpt.domain_om(1);
RAAN = (DynOpt.domain_RAAN(2)-DynOpt.domain_RAAN(1)).*rand(1) + DynOpt.domain_RAAN(1);
f0 = (DynOpt.domain_f0(2)-DynOpt.domain_f0(1)).*rand(1) + DynOpt.domain_f0(1);
T = (DynOpt.domain_T(2)-DynOpt.domain_T(1)).*rand(1) + DynOpt.domain_T(1);
DynOpt.orbit = [ecc; inclination; om; RAAN; f0; T];

% satellite init
[DynOpt, params,satellites_iner_ECI,satellites_attitude] = satellite_init_function(DynOpt, params, struct);
DynOpt.model_propagate = @model_propagate_local_v3;

%%% DEFINE MODEL NAME %%%
if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
    DynOpt.model = @InertialDynamicsIntegrator_V2_2_function;
elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_bias_v3_function;
elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_bias_v3_function;
    DynOpt.model_inertial = @InertialDynamicsIntegrator_V2_2_function;
end

% optimise input
DynOpt.optimise_input = struct.optimise_input;
DynOpt.optimise_params = struct.optimise_params;


%%%%%%%%%%%%%%%% LOAD FROM RL IF REINFORCEMENT LEARNING %%%%%%%%%%%%%%%%%%%
% set 0 Aw as default
DynOpt.Aw = zeros(3,1);

% SCALE FACTOR SETUP 
DynOpt.temp_scale = [1;1;1e-2;1e-3];
DynOpt.lambda = [1; 1; 1; 1];
DynOpt.y_weight = [1*ones(1,3) 1*ones(1,3*struct.nMagneto)];

%%%%% set orbit %%%%%
if DynOpt.generate_orbit
    DynOpt.orbit = struct.orbit;
end

if struct.RL
   DynOpt.RL = struct.RL_data;
   DynOpt.RL.S.satellites_attitude_singleopt = satellites_attitude;
   [DynOpt,satellites_iner_ECI,satellites_attitude] = RL_init_function_TD0(DynOpt,params,satellites_iner_ECI,struct);
else
    DynOpt.nMagneto = struct.nMagneto;
    DynOpt.dim_out = 3+3*DynOpt.nMagneto;
    DynOpt = scale_factor(DynOpt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLANT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DynOpt.generate_plant == 1
    % model simulation
    [DynOpt,params] = synthetic_integration_v3(DynOpt,params,struct,satellites_iner_ECI,satellites_attitude);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OBSERVER SETUP

if DynOpt.ObserverOn == 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% OPTIONS FOR GENERATE PLANT %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DynOpt.generate_plant == 1
        %%%%%%%%%%%%%%%%%%%%%%%% Derivatives setup %%%%%%%%%%%%%%%%%%%%%%%%
        if struct.print
            disp('Setting derivatives')
        end
        % first derivative 
        % number of derivatives
        DynOpt.c1_derivative = 1;
        %has to be smaller than DynOpt.Nts
        DynOpt.d1_derivative = 3;

        % buffer init
        DynOpt.buf_dY = zeros(9,DynOpt.d1_derivative);
        DynOpt.buf_dYhat = zeros(9,DynOpt.d1_derivative);

        % gradient buffer init
        DynOpt.buf_dyhat_grad = zeros(1,DynOpt.d1_derivative);
        
        % track error
        DynOpt.buf_trackerr = zeros(3,DynOpt.d1_derivative);

        % cost function derivative
        DynOpt.J_c1_derivative = 6;
        DynOpt.J_d1_derivative = 20;
        DynOpt.J_buf = zeros(1,DynOpt.J_d1_derivative);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%% Observer setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store state
        DynOpt.state = DynOpt.stateStory;  

        %true state and parameters
        DynOpt.Xtrue = [DynOpt.state(:,1); DynOpt.param_story(:,1)];
        DynOpt.aug_state_dim = length(DynOpt.Xtrue);

        % store quaternions in euler angles after simulation
        DynOpt.Opt_quat = zeros(3,1);

        %%%% init state - v1 %%%%
        offset = DynOpt.integration_pos*6;
        DynOpt.X_init = DynOpt.Xtrue;
        eul = [quat2eul(DynOpt.X_init(offset+1:offset+4)')';DynOpt.X_init(offset+5:end)];
        add_noise = 1*1e-3*rand(size(eul));
        eul_init = 1.3*eul+add_noise;
        quat_init = [eul2quat(eul_init(1:3)')'; eul_init(4:end)];
        % estimated attitude
        if DynOpt.noise_enable == 1
            DynOpt.X_init(offset+1:end) = quat_init;
        end
        DynOpt.X = DynOpt.X_init;
        DynOpt.Xtrue_init = DynOpt.Xtrue;

        % update params with the initial values
        DynOpt.params_update = @params_update_local_function;
        [params,DynOpt] = DynOpt.params_update(DynOpt.X_init,params,DynOpt);

        %%%%%%% storage %%%%%%%%
        DynOpt.WindowSamples = max(2,DynOpt.Nts*(DynOpt.w-1)+1);

        % measure derivative buffer
        DynOpt.Y =  zeros(9,DynOpt.w);
        DynOpt.Ytrue_full_story = [];
        DynOpt.Y_full_story = [];
        DynOpt.Yhat_full_story = [];

        % measure derivative buffer
        DynOpt.dY =  zeros(9,DynOpt.w);
        DynOpt.dY_full_story = [];
        DynOpt.dYhat_full_story = [];

        % measure integral buffer
        DynOpt.intY =  zeros(9,DynOpt.w);
        DynOpt.intY_full_story = [];
        DynOpt.intYhat_full_story = [];

        % buffer adaptive sampling
        DynOpt.Y_space = zeros(1,DynOpt.w);
        DynOpt.Y_space_full_story = 0;
        DynOpt.sensor_stop_thresh = 1.5;
        DynOpt.else_flag = 0;
        DynOpt.safety_density = struct.safety_density;
        DynOpt.safety_interval = int32(DynOpt.safety_density*DynOpt.WindowSamples);

        % dJ condition buffer (adaptive sampling)
        DynOpt.dJ_cond_story = [];

        % performance analysis
        DynOpt.performance_n_rows = 2*DynOpt.Nts;
        DynOpt.lambda_min = 0;
        DynOpt.measure_exp = 1;

        % corrupted trajectory
        DynOpt.Xstory(:,1) = DynOpt.X_init;

        % optimised trajectory
        DynOpt.OptXstory(:,1) = DynOpt.X;
        DynOpt.OptXstory_runtime(:,1) = DynOpt.X;

        % optimisation error
        DynOpt.OptErrorStory = DynOpt.OptXstory;

        % set reference trajectory
        if DynOpt.simulationModel == 1
            if ~isempty(DynOpt.param_estimate)
               DynOpt.OptXstoryTRUE = [DynOpt.state(1:DynOpt.StateDim,:);DynOpt.param_story];
            else
               DynOpt.OptXstoryTRUE = DynOpt.state(1:DynOpt.StateDim,:); 
            end
        else
            DynOpt.OptXstoryTRUE = [DynOpt.state(1,:);DynOpt.param_estimate'.*ones(2,length(DynOpt.time))];
        end

        % observer cost function init
        DynOpt.fcon_flag = struct.fcon_flag;
        DynOpt.J = 1e3;
        DynOpt.Jstory = DynOpt.J;
        DynOpt.Jdot_story = 0;
        DynOpt.J_meas = ones(9,1);
        DynOpt.J_der = ones(9,1);
        DynOpt.J_int = ones(9,1);
        DynOpt.J_quat = ones(9,1);
        DynOpt.temp_time = [];
        DynOpt.opt_chosen_time = [];
        DynOpt.grad_story = zeros(DynOpt.StateDim + length(DynOpt.param_estimate),1);

        % handle input
        DynOpt.U = zeros(3,length(DynOpt.time));
        DynOpt.recollect_input = 1;
        DynOpt.desired_attitude = zeros(3,length(DynOpt.time));
        DynOpt.e_flow = eye(length(DynOpt.X));
        params.DesiredAttitude = zeros(3,1);
        DynOpt.theta_u = [DynOpt.u_amp; DynOpt.d];

        % opt counters
        DynOpt.opt_counter = 0;
        DynOpt.select_counter = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        save plant_data
    else
        disp('ERROR LOADING PLANT DATA')
        disp('dont know how to manage workspaces in C')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT OBSERVER SETUP %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% OPTIMISATION SETUP %%%
    % constraints 
    DynOpt.Acon = [];
    DynOpt.Bcon = [];
    DynOpt.Acon_eq = [];
    DynOpt.Bcon_eq = [];
    DynOpt.nonlcon = @mycon;
    DynOpt.lb = [-Inf*ones(1,7), 0*ones(1,DynOpt.nparams)];
    DynOpt.ub = [Inf*ones(1,7), 1e-1*ones(1,DynOpt.nparams)];
    DynOpt.J_big = 0;
    % optimset
    DynOpt.TolX = 1e-10;
    DynOpt.TolFun = 1e-10;
    DynOpt.TolExit_J = struct.J_thresh;
    DynOpt.TolExit_X = 1e-10;
    DynOpt.outfun = @outputfcn;
    DynOpt.alphaVector = [1*ones(1,7), 1e-5*ones(1,DynOpt.nparams)];
    DynOpt.minAlgorithm = 'quasi-newton';
    DynOpt.GradObj = 'off';
    % define cost functions and measurements
    DynOpt.cost_function = @cost_function_v10;
    DynOpt.cost_function_name = 'cost_function_v10';
    if strcmp(struct.Observer,'EKF')
        DynOpt.get_measure = @get_measure_v6_function;
        DynOpt.get_measure_name = 'get_measure_v6';
    elseif strcmp(struct.Observer,'OPT')
        DynOpt.get_measure = @get_measure_v6_function;
        DynOpt.get_measure_name = 'get_measure_v6';
    end
    DynOpt.last_opt_time = 0;
    % max function counts
    DynOpt.maxFcount = struct.maxFcount;

    % set options
    if struct.print 
       DynOpt.display = 'iter';
    else
        DynOpt.display = 'off'; 
    end
    
    %%% for the gradient based functions %%%
    DynOpt.myoptioptions = optimset('Algorithm',DynOpt.minAlgorithm,'MaxIter', DynOpt.max_iter,'display',DynOpt.display,...
                                    'TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'MaxFunEvals',DynOpt.maxFcount,'OutputFcn',DynOpt.outfun,...
                                    'GradObj', DynOpt.GradObj);%...
%                                        'FinDiffType','central','FinDiffRelStep',DynOpt.alphaVector,'UseParallel',false);

    %%% for fminsearch only %%%
%     DynOpt.myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display',DynOpt.display,...
%                                     'TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'MaxFunEvals',DynOpt.maxFcount,'OutputFcn',DynOpt.outfun);


    DynOpt.fmin = struct.fmin;

    % gradient descent
    DynOpt.max_iter = struct.max_iter;
    DynOpt.opt_time = 0;

    % filtering options
    DynOpt.filter_flag = 0;
    DynOpt.filter_flag_opt = 0;
    DynOpt.filter_window = 10;
    
    % check update condition thresholds (if blue_flag==1 the previous are not checked)
    DynOpt.dJ_1 = struct.dJ_1;
    DynOpt.dJ_2 = struct.dJ_2;
    DynOpt.J_thresh = struct.J_thresh;
    DynOpt.Jdot_thresh = struct.Jdot_thresh;
    DynOpt.blue_flag = struct.blue_flag;
    
    % safety interval policy
    DynOpt.always_opt = struct.always_opt;
    DynOpt.flush_buffer = struct.flush_buffer;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % set RL flag
    DynOpt.RL_flag = 1;

    % memory load
    % save buffer mems
    if struct.load_mem
        
        start_step = 2;
        
        DynOpt.Y = struct.Y_last;
        DynOpt.dY = struct.dY_last;
        DynOpt.intY = struct.intY_last;
        DynOpt.buf_dY = struct.buf_dY_last;
        DynOpt.Y_full_story= struct.Y_full_story_last;
        DynOpt.dY_full_story = struct.dY_full_story_last;
        DynOpt.intY_full_story= struct.DynOpt.intY_full_story_last;

        DynOpt.buf_dYhat = struct.buf_dYhat_last;
        DynOpt.Yhat_full_story= struct.Yhat_full_story_last;
        DynOpt.dYhat_full_story = struct.dYhat_full_story_last;
        DynOpt.intYhat_full_story= struct.DynOpt.intYhat_full_story_last;

        DynOpt.Y_space = struct.Y_space_last;
        DynOpt.Y_space_full_story = struct.Y_space_last;
        
        DynOpt.OptXstoryTRUE = [struct.OptXstoryTRUE_last, DynOpt.OptXstoryTRUE(:,start_step:end)];
        DynOpt.OptXstory = [struct.OptXstory_last, DynOpt.OptXstory(:,2:end)];
        DynOpt.OptXstory_runtime = [struct.OptXstory_runtime_last, DynOpt.OptXstory_runtime(:,start_step:end)];
        DynOpt.Xstory = [struct.Xstory_last, DynOpt.Xstory(:,2:end)];
        
        DynOpt.time_tot = [struct.time_last, DynOpt.time(2:end)];
        
        DynOpt.input_true = [struct.input_true_last, DynOpt.input_true(:,2:end)];
        
        DynOpt.mag_field_story = [struct.mag_field_story_last, DynOpt.mag_field_story];
        
        DynOpt.eps_noise_story = [struct.eps_noise_story_last, DynOpt.eps_noise_story];
        
        DynOpt.buf_trackerr = struct.buf_trackerr_last;
        
        %%%% init state - v2 %%%%
        DynOpt.X_init = DynOpt.OptXstory_runtime(:,end);
        DynOpt.X = DynOpt.X_init;

        % update params with the initial values
        DynOpt.params_update = @params_update_local_function;
        [params,DynOpt] = DynOpt.params_update(DynOpt.X_init,params,DynOpt);
    end
    
%% Observer implementation 
    if DynOpt.OptimisationOn
        if strcmp(struct.Observer,'OPT')
            [DynOpt,params] = ObsOpt_bias_v6_function(DynOpt,params);
        else
            disp('Observer wrong option')
        end
    end
    
    
    %% save data
    save temp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error
    if DynOpt.simulationModel == 1
        DynOpt.OptErrorStory = DynOpt.OptXstoryTRUE - DynOpt.OptXstory;   
    end
    
    % measure
    DynOpt.Jstory = DynOpt.Jstory(:,2:end);
    DynOpt.Y_space_full_story = diff(DynOpt.Y_space_full_story);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Handle J chunks
    DynOpt.J_meas = DynOpt.J_meas(:,2:end);
    DynOpt.J_der = DynOpt.J_der(:,2:end);
    DynOpt.J_int = DynOpt.J_int(:,2:end);
    DynOpt.J_quat = DynOpt.J_quat(:,2:end);
    DynOpt.dJ_cond_story = DynOpt.dJ_cond_story(:,2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert quaterinions to euler angles 
    DynOpt.Opt_quat_runtime = DynOpt.wrap(quat2eul(DynOpt.OptXstory_runtime(DynOpt.integration_pos*6+1:DynOpt.integration_pos*6+4,DynOpt.past_length+1:end)'))';
    DynOpt.Opt_quat = DynOpt.wrap(quat2eul(DynOpt.OptXstory(DynOpt.integration_pos*6+1:DynOpt.integration_pos*6+4,DynOpt.past_length+1:end)'))';
    DynOpt.OptErrorStory_Euler = DynOpt.True_quat-DynOpt.Opt_quat_runtime;
    
    %%% tracking error %%%
    DynOpt.track_err = DynOpt.Opt_quat_runtime - DynOpt.target_attitude;
    for i = 1:size(DynOpt.track_err,2)
        for j = 1:size(DynOpt.track_err,1)
            [DynOpt.buf_trackerr, DynOpt.track_errdot(j,i)] = IterativePseudoDerivative(DynOpt.Ts,DynOpt.track_err(j,i),DynOpt.c1_derivative,DynOpt.d1_derivative,0,DynOpt.buf_trackerr);
        end
    end
      
    % performance
    DynOpt.lambda_min = DynOpt.lambda_min(2:end);
    DynOpt.grad_story = DynOpt.grad_story(:,2:end);
    
    % Optimisation time
    DynOpt.opt_time = DynOpt.opt_time(2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 1 
        
        %%% back step version 2
        back_step = 1;
        end_step = 1;
        
        %%% buffer length
        buf_len = size(DynOpt.buf_dY,2);
        story_len = size(DynOpt.Y_full_story,2);
        if story_len > buf_len
            DynOpt.buf_dY_last = DynOpt.Y_full_story(:,end-buf_len:end-1);
            DynOpt.buf_dYhat_last = DynOpt.Yhat_full_story(:,end-buf_len:end-1);
        else
            DynOpt.buf_dY_last = [zeros(size(DynOpt.Y_full_story,1),buf_len-story_len), DynOpt.Y_full_story];
            DynOpt.buf_dYhat_last = [zeros(size(DynOpt.Y_full_story,1),buf_len-story_len), DynOpt.Yhat_full_story];
        end
        
        DynOpt.Y_last = DynOpt.Y;
        DynOpt.dY_last = DynOpt.dY;
        DynOpt.intY_last = DynOpt.intY;
        
        DynOpt.Y_full_story_last = DynOpt.Y_full_story(:,back_step:end-end_step);
        DynOpt.dY_full_story_last = DynOpt.dY_full_story(:,back_step:end-end_step);
        DynOpt.intY_full_story_last = DynOpt.intY_full_story(:,back_step:end-end_step);

        
        DynOpt.Yhat_full_story_last = DynOpt.Yhat_full_story(:,back_step:end-end_step);
        DynOpt.dYhat_full_story_last = DynOpt.dYhat_full_story(:,back_step:end-end_step);
        DynOpt.intYhat_full_story_last = DynOpt.intYhat_full_story(:,back_step:end-end_step);
        
        DynOpt.Y_space_last = DynOpt.Y_space;
        DynOpt.Y_space_full_story_last = DynOpt.Y_space_full_story;
        
        DynOpt.OptXstoryTRUE_last = DynOpt.OptXstoryTRUE(:,back_step:end);
        DynOpt.OptXstory_last = DynOpt.OptXstory(:,back_step:end);
        DynOpt.OptXstory_runtime_last = DynOpt.OptXstory_runtime(:,back_step:end);
        DynOpt.Xstory_last = DynOpt.Xstory(:,back_step:end);
        
        DynOpt.time_last = DynOpt.time_tot;
        
        % keep end because it starts with zeros
        DynOpt.input_true_last = DynOpt.input_true(:,back_step:end);
        
        DynOpt.mag_field_story_last = DynOpt.mag_field_story(:,back_step:end-end_step);
        
        DynOpt.eps_noise_story_last = DynOpt.eps_noise_story(:,back_step:end);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

