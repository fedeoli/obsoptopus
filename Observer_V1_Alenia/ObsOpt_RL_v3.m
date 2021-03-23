function [DynOpt, params] = ObsOpt_RL_v3(struct,RL)

%% Init Section
close all
% clc

% dependencies
addpath(genpath([pwd '/Lib']));
addpath(genpath([pwd '/Satellite']));

% new random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% global data structure
global DynOpt params

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

%%% init flags and setup %%%
% model name
DynOpt.modelname = struct.model;

% random init
DynOpt.rand_init = struct.rand_init;

% observer setup
DynOpt.ObserverOn = struct.ObserverOn;
DynOpt.OptimisationOn = struct.OptimisationOn;
DynOpt.scale_factor = struct.scale_factor;

% starts finding the state/parameters backward, -1, (actual time t_k) or Forward, 1 (i.e. at time  t_{k-(w-1)*Nts}
DynOpt.ForwardOptimization = struct.forward; 

% measure noise
DynOpt.noise_enable = struct.noise_enable;
DynOpt.measure_amp = struct.noise_amp;
DynOpt.nMagneto = struct.nMagneto;

% simulate the plant
DynOpt.simulationModel = struct.simulationModel;

% wrap function
DynOpt.wrap = @unwrap;

% print option
DynOpt.print = struct.print;

% montecarlo flag
DynOpt.montecarlo = struct.montecarlo;
if DynOpt.montecarlo == 1
    % state bounds
    DynOpt.params_init = struct.params_init;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLANT model and synthetic data

%%%%%%%%%%%%%%%%%%% SET INITIAL CONDITIONS FOR GREEDY %%%%%%%%%%%%%%%%%
DynOpt.u_freq = RL.S.F0(RL.S.i);
DynOpt.target_attitude = RL.S.T0(:,RL.S.i);
DynOpt.u_amp = RL.S.A(1,RL.S.i);
DynOpt.d = RL.S.A(2,RL.S.i);

%%% INPUT SETUP%%%
DynOpt.set_input = @set_input_v5;
DynOpt.switch_pwm = 1;
DynOpt.t_lowpass = 1;
DynOpt.lowpass_pwm = struct.lowpass_pwm;

% satellite init
satellite_init;
DynOpt.model_propagate = @model_propagate_local_v3;

%%% DEFINE MODEL NAME %%%
if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
    DynOpt.model = @InertialDynamicsIntegrator_V2_2;
elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_bias_v2;
    DynOpt.model_flow = @satellite_dflow_input_v1;
elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_bias_v2;
    DynOpt.model_inertial = @InertialDynamicsIntegrator_V2_2;
end

%%% initial condition %%%
params.SatellitesAttitude = [eul2quat(RL.S.S0(1:3,RL.S.i)')';RL.S.S0(4:6,RL.S.i)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLANT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DynOpt.generate_plant == 1
    % model simulation
    tic
    synthetic_integration_v2;
    DynOpt.synthetic_time = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OBSERVER SETUP

if DynOpt.ObserverOn == 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% OPTIONS FOR GENERATE PLANT %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DynOpt.generate_plant == 1
        %%%%%%%%%%%%%%%%%%%%%%%% Derivatives setup %%%%%%%%%%%%%%%%%%%%%%%%
        if DynOpt.print
            disp('Setting derivatives')
        end
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%% Observer setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store state
        DynOpt.state = DynOpt.stateStory;  

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
        offset = DynOpt.integration_pos*6;
        if DynOpt.noise_enable == 1
            % initial error on EA instead of quaternion
            if struct.ErrorOnEA == 1
                if DynOpt.integration_pos == 1
                    X_init = [DynOpt.Xtrue(1:offset); quat2eul(DynOpt.Xtrue(offset+1:offset+4,:)')'; DynOpt.Xtrue(offset+5:end,:)];
                else
                    X_init = [quat2eul(DynOpt.Xtrue(offset+1:offset+4,:)')'; DynOpt.Xtrue(offset+5:end,:)];
                end
                init_error_amp = DynOpt.init_error_amp(1:end);
                stateDim = DynOpt.StateDim-1;
            else
                X_init = DynOpt.Xtrue;
                init_error_amp = DynOpt.init_error_amp;
                stateDim = DynOpt.StateDim;
            end
            if DynOpt.rand_init == 1
                if struct.ErrorOnEA == 0
                    state_init = (1+1e-2*init_error_amp.*randn(stateDim,1)).*X_init(1:stateDim);
                else
                    state_init = [X_init(1:3)+init_error_amp(1:3).*randn(3,1);(1+1e-2*init_error_amp(4:end).*randn(stateDim-3,1)).*X_init(4:stateDim)];
                end
                param_init = (1+1e-2*DynOpt.init_param_error_amp.*randn(length(DynOpt.param_estimate),1)).*DynOpt.Xtrue(DynOpt.StateDim+1:DynOpt.StateDim+DynOpt.nparams);
            else
                if struct.ErrorOnEA == 0
                    state_init = (1+1e-2*init_error_amp.*ones(stateDim,1)).*X_init(1:stateDim);
                else
                    init_error_amp(DynOpt.integration_pos*6+1:DynOpt.integration_pos*6+3) = deg2rad(init_error_amp(DynOpt.integration_pos*6+1:DynOpt.integration_pos*6+3));
                    state_init = [X_init(offset+1:offset+3)+init_error_amp(offset+1:offset+3).*ones(3,1);(1+1e-2*init_error_amp(offset+5:end).*ones(stateDim-(3+offset),1)).*X_init(offset+4:stateDim)];
                    if DynOpt.integration_pos == 1
                        pos_init = (1+1e-2*init_error_amp(1:6).*ones(6,1)).*X_init(1:6);
                        state_init = [pos_init; state_init];
                    end
                end
                
                param_init = (1+1e-2*transpose(DynOpt.init_param_error_amp).*ones(length(DynOpt.param_estimate),1)).*DynOpt.Xtrue(DynOpt.StateDim+1:DynOpt.StateDim+DynOpt.nparams);
            end
            if DynOpt.paramsAllPos
                param_init = abs(param_init);
            end
            if DynOpt.input_tuning == 1
               param_init = [param_init; DynOpt.param_story(end-2:end,1)]; 
            end
            if struct.ErrorOnEA == 1
                init_quat = eul2quat(state_init(offset+1:offset+3)')';
                if DynOpt.integration_pos == 1
                    state_init = [state_init(1:6); init_quat; state_init(offset+4:end)];
                else
                    state_init = [init_quat; state_init(4:end)];
                end
                
            end
            DynOpt.X  = [state_init; param_init];   
        else
            DynOpt.X = DynOpt.Xtrue; 
        end

        DynOpt.X_init = DynOpt.X;
        DynOpt.Xtrue_init = DynOpt.Xtrue;

        % update params with the initial values
        params_update(DynOpt.X_init);

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
        DynOpt.sensor_stop_thresh = 1.5;
        DynOpt.theta = struct.theta;
        DynOpt.beta = struct.beta;
        DynOpt.gamma = struct.gamma;
        DynOpt.else_flag = 0;
        DynOpt.safety_density = struct.safety_density;
        DynOpt.safety_interval = int32(DynOpt.safety_density*DynOpt.WindowSamples);

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
        DynOpt.OptXstory_runtime = zeros(length(DynOpt.X),length(DynOpt.time));
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
        DynOpt.J = 1e3;
        DynOpt.Jstory = DynOpt.J;
        DynOpt.Jdot_story = 0;
        DynOpt.J_meas = ones(DynOpt.dim_out,1);
        DynOpt.J_der = ones(DynOpt.dim_out,1);
        DynOpt.J_int = ones(DynOpt.dim_out,1);
        DynOpt.J_dyn = ones(DynOpt.dim_out,1);
        DynOpt.J_quat = ones(DynOpt.dim_out,1);
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
        try
            save temp
            load plant_data
            load temp
            delete temp.mat
        catch
            disp('ERROR LOADING PLANT DATA')
            return
        end
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
    DynOpt.cost_function = @cost_function_v7;
    DynOpt.cost_function_name = 'cost_function_v7';
    if strcmp(struct.Observer,'EKF')
        DynOpt.get_measure = @get_measure_v5;
        DynOpt.get_measure_name = 'get_measure_v5';
    elseif strcmp(struct.Observer,'OPT')
        DynOpt.get_measure = @get_measure_v5;
        DynOpt.get_measure_name = 'get_measure_v5';
    end
    DynOpt.params_update = @params_update_local;
    DynOpt.last_opt_time = 0;
    % multistart
    DynOpt.multistart = struct.multistart;
    DynOpt.globalsearch = struct.globalsearch;
    if DynOpt.multistart == 1
        DynOpt.nstart = struct.nstart;
        DynOpt.fmin_name = 'fmincon';
    end
    
    if DynOpt.print == 1
       DynOpt.display = 'iter';
    else
       DynOpt.display = 'off';
    end

    % set options
    if DynOpt.fcon_flag == 0
        DynOpt.myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display',DynOpt.display,'TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
%         myoptioptions = optimset('Algorithm','quasi-newton','MaxIter', DynOpt.max_iter,'display','iter','GradConstr','off','GradObj','off','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        DynOpt.fmin = struct.fmin;
    else
%         myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display','iter','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        DynOpt.myoptioptions = optimset('Algorithm','interior-point','MaxIter', DynOpt.max_iter,'display',DynOpt.display,'GradConstr','off','GradObj','off','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
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
    
    % check update condition thresholds (if blue_flag==1 the previous are not checked)
    DynOpt.dJcond_thresh = struct.dJcond_thresh;
    DynOpt.J_thresh = struct.J_thresh;
    DynOpt.Jdot_thresh = struct.Jdot_thresh;
    DynOpt.blue_flag = struct.blue_flag;
    
    % optimise input
    DynOpt.optimise_input = struct.optimise_input;
    DynOpt.optimise_params = struct.optimise_params;
    
    % safety interval policy
    DynOpt.always_opt = struct.always_opt;
    DynOpt.flush_buffer = struct.flush_buffer;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Observer implementation 
    if DynOpt.OptimisationOn
        if strcmp(struct.Observer,'OPT')
            ObsOpt_bias_v5;
        elseif strcmp(struct.Observer,'EKF')
            ObsOpt_EKF_v1;
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
    DynOpt.Y_full_story = DynOpt.Y_full_story(:,2:end);
    DynOpt.Yhat_full_story = DynOpt.Yhat_full_story(:,2:end);
    DynOpt.dY_full_story = DynOpt.dY_full_story(:,2:end);
    DynOpt.dYhat_full_story = DynOpt.dYhat_full_story(:,2:end);
    DynOpt.Jstory = DynOpt.Jstory(:,2:end);
    DynOpt.Y_space_full_story = diff(DynOpt.Y_space_full_story);
    DynOpt.intY_full_story = DynOpt.intY_full_story(:,2:end);
    DynOpt.intYhat_full_story = DynOpt.intYhat_full_story(:,2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Handle J chunks
    DynOpt.J_meas = DynOpt.J_meas(:,2:end);
    DynOpt.J_der = DynOpt.J_der(:,2:end);
    DynOpt.J_int = DynOpt.J_int(:,2:end);
    DynOpt.J_dyn = DynOpt.J_dyn(:,2:end);
    DynOpt.J_quat = DynOpt.J_quat(:,2:end);
    DynOpt.dJ_cond_story = DynOpt.dJ_cond_story(:,2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert quaterinions to euler angles 
    DynOpt.Opt_quat = DynOpt.wrap(quat2eul(DynOpt.OptXstory(DynOpt.integration_pos*6+1:DynOpt.integration_pos*6+4,:)'))';
    DynOpt.OptErrorStory_Euler = DynOpt.True_quat-DynOpt.Opt_quat;
      
    % performance
    DynOpt.lambda_min = DynOpt.lambda_min(2:end);
%     DynOpt.Jstory = DynOpt.Jstory(1,2:end);
    DynOpt.grad_story = DynOpt.grad_story(:,2:end);
    
    % Optimisation time
    DynOpt.opt_time = DynOpt.opt_time(2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

