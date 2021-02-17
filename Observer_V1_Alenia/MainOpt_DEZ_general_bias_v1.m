function [DynOpt, params] = MainOpt_DEZ_general_bias_v1(struct)

%% Init Section
close all
clc

% dependencies
addpath(genpath([pwd '/Lib']));
addpath(genpath([pwd '/Satellite']));

% new random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% global data structure
global DynOpt params

% identify
DynOpt.identify = struct.identify;

%%% init flags and setup %%%
% model name
DynOpt.modelname = struct.model;

% observer setup
DynOpt.ObserverOn = struct.ObserverOn;
DynOpt.OptimisationOn = struct.OptimisationOn;
DynOpt.scale_factor = struct.scale_factor;

% starts finding the state/parameters backward, -1, (actual time t_k) or Forward, 1 (i.e. at time  t_{k-(w-1)*Nts}
DynOpt.ForwardOptimization = struct.forward; 

%%% ENABLE INPUT TUNING MODE %%%
DynOpt.input_tuning = struct.input_tuning;

% Inter-sampling Time, number of samples within intersampled data
DynOpt.Nts = struct.Nts; 

%number of inter-sampled data in a window: the total number of sampled data in a window are (DynOpt.w-1)*DynOpt.Nts+1
DynOpt.w = struct.w;

% measure noise
DynOpt.noise_enable = struct.noise_enable;
DynOpt.measure_amp = struct.noise_amp;

% simulate the plant
DynOpt.simulationModel = struct.simulationModel;

% montecarlo flag
DynOpt.montecarlo = struct.montecarlo;
if DynOpt.montecarlo == 1
    % state bounds
    DynOpt.params_init = struct.params_init;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLANT model and synthetic data
% satellite init
satellite_init;

%%% INPUT SETUP%%%
DynOpt.u_amp = struct.u_amp;
DynOpt.u_freq = struct.u_freq;
DynOpt.d = struct.d;
DynOpt.target_attitude = struct.target_attitude;
DynOpt.set_input = @set_input_v4;
DynOpt.switch_pwm = 1;
DynOpt.t_lowpass = 1;

%%% DEFINE MODEL NAME %%%
if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
    DynOpt.model = @InertialDynamicsIntegrator_V2_2;
elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_bias_v2;
    DynOpt.model_flow = @satellite_dflow_input;
elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_bias_v2;
    DynOpt.model_inertial = @InertialDynamicsIntegrator_V2_2;
end

% bias flag
DynOpt.bias_model = 1;

% model simulation
tic
synthetic_integration;
toc

%% Derivatives setup

disp('Setting derivatives')
% first derivative 
% number of derivatives
DynOpt.c1_derivative = 1;
%has to be smaller than DynOpt.Nts
DynOpt.d1_derivative = 3;

% buffer init
DynOpt.buf_dy = zeros(DynOpt.dim_out,DynOpt.d1_derivative);
DynOpt.buf_dyhat = zeros(DynOpt.dim_out,DynOpt.d1_derivative);

% gradient buffer init
DynOpt.buf_dyhat_grad = zeros(1,DynOpt.d1_derivative);

% cost function derivative
DynOpt.J_c1_derivative = 6;
DynOpt.J_d1_derivative = 20;
DynOpt.J_buf = zeros(1,DynOpt.J_d1_derivative);

% store state
DynOpt.state = DynOpt.stateStory;


%% OBSERVER SETUP

if DynOpt.ObserverOn == 1
    %true state and parameters
    DynOpt.Xtrue = [DynOpt.state(:,1); DynOpt.param_story(:,1)];
    DynOpt.aug_state_dim = length(DynOpt.Xtrue);
    
    % store quaternions in euler angles after simulation
    DynOpt.Opt_quat = zeros(3,1);
    
    % init state
    DynOpt.init_error_amp = struct.init_error_amp;
    DynOpt.noise_enable = struct.noise_enable;
    DynOpt.init_param_error_amp = struct.init_param_error_amp;
    DynOpt.paramsAllPos = struct.ParamsAllPos;
    if DynOpt.noise_enable == 1
        state_init = (1+1e-2*DynOpt.init_error_amp.*randn(DynOpt.StateDim,1)).*DynOpt.Xtrue(1:DynOpt.StateDim);
        param_init = (1+1e-2*DynOpt.init_param_error_amp.*randn(length(DynOpt.param_estimate),1)).*DynOpt.Xtrue(DynOpt.StateDim+1:DynOpt.StateDim+DynOpt.nparams);
        if DynOpt.paramsAllPos
            param_init = abs(param_init);
        end
        if DynOpt.input_tuning == 1
           param_init = [param_init; DynOpt.param_story(end-2:end,1)]; 
        end
        DynOpt.X  = [state_init; param_init];   
    else
        DynOpt.X = DynOpt.Xtrue; 
    end
    
    DynOpt.X_init = DynOpt.X;
    DynOpt.Xtrue_init = DynOpt.Xtrue;
    
    % update params with the initial values
    params_update(DynOpt.X_init);
    
    % observer optimisation options
    DynOpt.fcon_flag = struct.fcon_flag;
    DynOpt.max_iter = struct.max_iter;
    DynOpt.identify = struct.identify;
    
    %%%% OPTIMISATION SETUP %%%
    % constraints 
    DynOpt.Acon = [];
    DynOpt.Bcon = [];
    DynOpt.Acon_eq = [];
    DynOpt.Bcon_eq = [];
    DynOpt.lb = struct.lb;
    DynOpt.ub = struct.ub;
    DynOpt.nonlcon = [];
    DynOpt.J_big = 0;
    % optimset
    DynOpt.TolX = 1e-10;
    DynOpt.TolFun = 1e-10;
    DynOpt.TolExit_J = struct.J_thresh;
    DynOpt.TolExit_X = 1e-10;
    DynOpt.outfun = @outputfcn;
    % define cost functions and measurements
    if DynOpt.input_tuning == 0
        DynOpt.cost_function = @cost_function_v4;
        DynOpt.cost_function_name = 'cost_function_v4';
    else
        DynOpt.cost_function = @cost_function_input_v1;
        DynOpt.cost_function_name = 'cost_function_input_v1';
    end
    DynOpt.get_measure = @get_measure_v3;
    DynOpt.get_measure_name = 'get_measure_v3';
    DynOpt.model_propagate = @model_propagate_local;
    DynOpt.params_update = @params_update_local;
    DynOpt.last_opt_time = 0;
    % multistart
    DynOpt.multistart = struct.multistart;
    DynOpt.globalsearch = struct.globalsearch;
    if DynOpt.multistart == 1
        DynOpt.nstart = struct.nstart;
        DynOpt.fmin_name = 'fmincon';
    end

    % set options
    if DynOpt.fcon_flag == 0
        DynOpt.myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display','iter','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
%         myoptioptions = optimset('Algorithm','quasi-newton','MaxIter', DynOpt.max_iter,'display','iter','GradConstr','off','GradObj','off','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        DynOpt.fmin = struct.fmin;
    else
%         myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display','iter','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        DynOpt.myoptioptions = optimset('Algorithm','interior-point','MaxIter', DynOpt.max_iter,'display','iter','GradConstr','off','GradObj','off','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        DynOpt.fmin = struct.fmin;
    end
  
    % gradient descent
    DynOpt.y_end = struct.y_end;
    DynOpt.alpha_grad = struct.alpha_grad;
    DynOpt.max_iter = struct.max_iter;
    DynOpt.grad_thresh = struct.grad_thresh;
    DynOpt.opt_time = 0;

    % filtering options
    DynOpt.filter_flag = 0;
    DynOpt.filter_flag_opt = 0;
    DynOpt.filter_window = 10;

    %%%%%%% storage %%%%%%%%
    DynOpt.WindowSamples = max(2,DynOpt.Nts*(DynOpt.w-1)+1);
    
    % measure derivative buffer
    DynOpt.Y =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.Ytrue_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.Y_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.Yhat_full_story = zeros(DynOpt.dim_out,1);
    
    % measure derivative buffer
    DynOpt.dY =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.dY_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.dYhat_full_story = zeros(DynOpt.dim_out,1);
    
    % measure integral buffer
    DynOpt.intY =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.intY_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.intYhat_full_story = zeros(DynOpt.dim_out,1);
    
    % buffer adaptive sampling
    DynOpt.Y_space = zeros(1,DynOpt.w);
    DynOpt.Y_space_full_story = 0;
    DynOpt.dJcond_thresh = struct.dJcond_thresh;
    DynOpt.sensor_stop_thresh = 1.5;
    DynOpt.theta = struct.theta;
    DynOpt.beta = struct.beta;
    DynOpt.gamma = struct.gamma;
    DynOpt.else_flag = 0;
    DynOpt.safety_density = struct.safety_density;
    DynOpt.safety_interval = int32(DynOpt.safety_density*DynOpt.WindowSamples);
    
    %%% IF INPUT TUNING MODE - FIXED SAMPLING %%%
    if DynOpt.input_tuning == 1
        DynOpt.dJcond_thresh = 0;
    end
    
    % dJ condition buffer (adaptive sampling)
    DynOpt.dJ_cond_story = zeros(4,1);
    
    % performance analysis
    DynOpt.performance_n_rows = 2*DynOpt.Nts;
    DynOpt.lambda_min = 0;
    DynOpt.measure_exp = 1;

    % corrupted trajectory
    DynOpt.Xstory = zeros(length(DynOpt.X),length(DynOpt.time));
    DynOpt.Xstory(:,1) = DynOpt.X_init;

    % optimised trajectory
    DynOpt.OptXstory = zeros(length(DynOpt.X),length(DynOpt.time));
    DynOpt.OptXstory(:,1) = DynOpt.X;

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
    DynOpt.J = 1e3;
    DynOpt.Jstory = DynOpt.J;
    DynOpt.Jdot_story = 0;
    DynOpt.J_meas = 1;
    DynOpt.J_der = 1;
    DynOpt.J_int = 1;
    DynOpt.temp_time = [];
    DynOpt.opt_chosen_time = [];
    DynOpt.grad_story = zeros(DynOpt.StateDim + length(DynOpt.param_estimate),1);


    % check update condition thresholds (if blue_flag==1 the previous are not checked)
    DynOpt.J_thresh = struct.J_thresh;
    DynOpt.Jdot_thresh = struct.Jdot_thresh;
    DynOpt.blue_flag = struct.blue_flag;
    
    % handle input
    DynOpt.U = zeros(3,length(DynOpt.time));
    DynOpt.recollect_input = 1;
    DynOpt.optimise_input = struct.optimise_input;
    DynOpt.optimise_params = struct.optimise_params;
    DynOpt.desired_attitude = zeros(3,length(DynOpt.time));
    DynOpt.e_flow = eye(length(DynOpt.X));
    params.DesiredAttitude = zeros(3,1);
    DynOpt.theta_u = [DynOpt.u_amp; DynOpt.d];
    
    % safety interval policy
    DynOpt.always_opt = struct.always_opt;
    DynOpt.flush_buffer = struct.flush_buffer;
    
    % opt counters
    DynOpt.opt_counter = 0;
    DynOpt.select_counter = 0;
    
%% Observer implementation 
    if DynOpt.OptimisationOn
        ObsOpt_bias_v3;
    end
    
    
    %% save data
    save temp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error
    if DynOpt.simulationModel == 1
        DynOpt.OptErrorStory = DynOpt.OptXstoryTRUE - DynOpt.OptXstory;   
    end
    
    % measure
    DynOpt.Y_full_story = DynOpt.Y_full_story(:,2:end);
    DynOpt.Yhat_full_story = DynOpt.Yhat_full_story(:,2:end);
    DynOpt.dY_full_story = DynOpt.dY_full_story(:,2:end);
    DynOpt.dYhat_full_story = DynOpt.dYhat_full_story(:,2:end);
    DynOpt.Jstory = DynOpt.Jstory(:,2:end);
    DynOpt.Y_space_full_story = diff(DynOpt.Y_space_full_story);
    DynOpt.intY_full_story = DynOpt.intY_full_story(:,2:end);
    DynOpt.intYhat_full_story = DynOpt.intYhat_full_story(:,2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert quaterinions to euler angles 
    DynOpt.Opt_quat = unwrap(RotationConversion_V2_1('QtoEA321', DynOpt.OptXstory(1:4,:)')*pi/180)';
      
    % performance
    DynOpt.lambda_min = DynOpt.lambda_min(2:end);
    DynOpt.Jstory = DynOpt.Jstory(1,2:end);
    DynOpt.grad_story = DynOpt.grad_story(:,2:end);
    
    % Optimisation time
    DynOpt.opt_time = DynOpt.opt_time(2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

