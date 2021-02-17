function [DynOpt, params] = MainOpt_DEZ_general_v22_fun_params_alenia(struct)

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

%% PLANT model and data
% satellite init
satellite_init;

%%% INPUT SETUP%%%
DynOpt.u_amp = struct.u_amp;
DynOpt.u_freq = struct.u_freq;
DynOpt.d = struct.d;
DynOpt.target_attitude = struct.target_attitude;
DynOpt.set_input = @set_input_v4;

%%% DEFINE MODEL NAME %%%
if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
    DynOpt.model = @InertialDynamicsIntegrator_V2_2;
elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
    DynOpt.model = @AttitudeDynamics_V2_2;
    DynOpt.model_flow = @satellite_dflow_input;
elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
    DynOpt.model = [@InertialDynamicsIntegrator_V2_2, @AttitudeDynamics_V2_2];
end

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
        DynOpt.cost_function_name = 'cost_function_v2';
    else
        DynOpt.cost_function = @cost_function_input_v1;
        DynOpt.cost_function_name = 'cost_function_input_v1';
    end
    DynOpt.get_measure = @get_measure_v3;
    DynOpt.get_measure_name = 'get_measure_v3';
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
        myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display','iter','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
%         myoptioptions = optimset('Algorithm','quasi-newton','MaxIter', DynOpt.max_iter,'display','iter','GradConstr','off','GradObj','off','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        DynOpt.fmin = struct.fmin;
    else
%         myoptioptions = optimset('MaxIter', DynOpt.max_iter,'display','iter','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
        myoptioptions = optimset('Algorithm','interior-point','MaxIter', DynOpt.max_iter,'display','iter','GradConstr','off','GradObj','off','TolFun',DynOpt.TolFun,'TolX',DynOpt.TolX,'OutputFcn',DynOpt.outfun);
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
    DynOpt.desired_attitude = zeros(3,length(DynOpt.time));
    DynOpt.e_flow = eye(length(DynOpt.X));
    params.DesiredAttitude = zeros(3,1);
    DynOpt.theta_u = [DynOpt.u_amp; DynOpt.d];
    
    % safety interval policy
    DynOpt.always_opt = struct.always_opt;
    
    % opt counters
    DynOpt.opt_counter = 0;
    DynOpt.select_counter = 0;
    
%% Observer implementation 
    if DynOpt.OptimisationOn
        disp('Processing data with the optimization-based observer...')
        run_time = tic;
        for k=1:length(DynOpt.time)

            % update actual index
            DynOpt.ActualTimeIndex = k;

            % reference state - used for noise
            DynOpt.Xtrue = [DynOpt.state(:,k);DynOpt.param_story(:,k)];
            
            %forward propagation of the previous estimate
            if(k>1)
                % INTEGRATION OF BOTH POSITION AND ATTITUDE - STACK 
                if DynOpt.input_tuning == 0
                    % Control allocation inside "params" structure          
                    [DynOpt.X, params] = model_propagate_local(DynOpt.ActualTimeIndex,DynOpt.Ts,DynOpt.OptXstory(:,DynOpt.ActualTimeIndex-1),params);
                    DynOpt.OptXstory(:,k) = DynOpt.X; 
                else
                    % current time index because it uses the true state
                    temp_state = DynOpt.OptXstoryTRUE(1:DynOpt.StateDim+DynOpt.nparams,k);
                    DynOpt.X = [temp_state; DynOpt.X(end-2:end)];
                    DynOpt.OptXstory(:,k) = DynOpt.X; 
                end
            end

            %%%%%%%%% MEASUREMENT %%%%%%%%%%%
            % read measure 
            measure_forward = 1;

            if DynOpt.simulationModel == 1
                [DynOpt.buf_dy,Y_true] = DynOpt.get_measure(DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dy,DynOpt.intY_full_story,params);
                DynOpt.measure_noise = DynOpt.measure_amp*randn(length(params.observed_state),1);
                % copy to Y noise and corrupt only the measure 
                Y_noise = Y_true;
                Y_noise(:,1) = Y_noise(:,1) + [DynOpt.measure_noise.*ones(length(params.observed_state),1); zeros(DynOpt.dim_out-length(params.observed_state),1)];
            else
                Y_true = EvaluateCostFunctionOnWindow_Output_general_data(k,measure_forward);
                Y_noise = Y_true;
            end

            % no filtering
            Y_filter = Y_noise;

            % store total memory
            DynOpt.Ytrue_full_story(:,end+1) = Y_true(:,1);
            DynOpt.Y_full_story(:,end+1) = Y_filter(:,1);
            DynOpt.dY_full_story(:,end+1) = Y_filter(:,2);
            DynOpt.intY_full_story(:,end+1) = Y_filter(:,3);
            
            % fisrt bunch of data - read Y every Nts and check if the signal is
            dJ_cond(DynOpt.theta,DynOpt.beta,DynOpt.gamma);
            distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);
            DynOpt.distance_safe_flag = (distance < DynOpt.safety_interval);
            % SENSOR SATURATION/ACCURACY
            DynOpt.sensor_stop_flag = (DynOpt.Y_full_story(DynOpt.ActualTimeIndex) < DynOpt.sensor_stop_thresh);
            % if last sample is too near or the condition isn't met then
            % propagate. However, check that the last optimisation has been
            % performed at most WindowSamples times ago
            if  ((distance < DynOpt.Nts) || (DynOpt.dJ_cond < DynOpt.dJcond_thresh)) && DynOpt.distance_safe_flag
                %%%% ESTIMATED measurements
                % measures
                [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params);
                DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
                DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);

                % clean 
                clc
            else

                % Display iteration slengthtep
                if DynOpt.input_tuning == 0
                    disp(['TARGET: I1: ', num2str(DynOpt.OptXstoryTRUE(end-2,DynOpt.ActualTimeIndex)), ' I2: ', num2str(DynOpt.OptXstoryTRUE(end-1,DynOpt.ActualTimeIndex)), ' I3: ', num2str(DynOpt.OptXstoryTRUE(end,DynOpt.ActualTimeIndex))])
                    disp(['INIT: I1: ', num2str(DynOpt.X_init(8)), ' I2: ', num2str(DynOpt.X_init(9)), ' I3: ', num2str(DynOpt.X_init(10))])
                    disp(['CURRENT: I1: ', num2str(params.sat(1).I(1,1)), ' I2: ', num2str(params.sat(1).I(2,2)), ' I3: ', num2str(params.sat(1).I(3,3))])
                else
                    disp(['INIT: Au: ', num2str(DynOpt.X_init(11)), ' d: ', num2str(DynOpt.X_init(12)), ' f: ', num2str(DynOpt.X_init(13))])
                    disp(['CURRENT: Au: ', num2str(DynOpt.X(11)), ' d: ', num2str(DynOpt.X(12)), ' f: ', num2str(DynOpt.X(13))])
                end
                disp(['n window: ', num2str(DynOpt.w),'  n samples: ', num2str(DynOpt.Nts)])
                disp(['Iteration Number: ', num2str(k),'/',num2str(length(DynOpt.time))])
                disp(['Last cost function: ', num2str(DynOpt.Jstory(end))]);
                disp(['N. optimisations RUN: ',num2str(DynOpt.opt_counter)]);
                disp(['N. optimisations SELECTED: ',num2str(DynOpt.select_counter)]);

                %%%% OUTPUT measurements - buffer of w elements
                % measures
                DynOpt.Y(:,1:end-1) = DynOpt.Y(:,2:end);
                DynOpt.Y(:,end) = Y_filter(:,1);

                % measures derivative
                DynOpt.dY(:,1:end-1) = DynOpt.dY(:,2:end);
                DynOpt.dY(:,end) = Y_filter(:,2);
                
                % measures integral
                DynOpt.intY(:,1:end-1) = DynOpt.intY(:,2:end);
                DynOpt.intY(:,end) = Y_filter(:,3);
                
                % backup
                Y_space_backup = DynOpt.Y_space;
                Y_space_full_story_backup = DynOpt.Y_space_full_story;
                
                % adaptive sampling
                DynOpt.Y_space(1:end-1) = DynOpt.Y_space(2:end);
                DynOpt.Y_space(end) = DynOpt.ActualTimeIndex;
                DynOpt.Y_space_full_story(end+1) = DynOpt.ActualTimeIndex;

                % store measure times
                DynOpt.temp_time = [DynOpt.temp_time k];

                if (k < max(1,DynOpt.WindowSamples)) 
                    [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params);
                    DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                    DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
                    DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);
                else    
                    
                    % measures
                    [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params);
                    DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                    DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
                    DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    %%% forward optimization %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if(DynOpt.ForwardOptimization == 1) 

                        % back time index
                        n_samples = min(length(DynOpt.Y_space_full_story)-1,DynOpt.w);
                        buf_dist = diff(DynOpt.Y_space_full_story(end-n_samples:end));
                        DynOpt.BackTimeIndex = max(1,k-sum(buf_dist)); 

                        % set of initial conditions
                        DynOpt.temp_x0 = DynOpt.OptXstory(:,DynOpt.BackTimeIndex);
                        
                        %%% INPUT TUNING - HYBRID INIT %%%
                        if DynOpt.input_tuning == 1
                            input_state = DynOpt.temp_x0(end-2:end-1);
                        end
                        
                        % Optimisation
                        % save current state for recovery
                        recovery_pos = params.SatellitesCoordinates;
                        recovery_att = params.SatellitesAttitude;

                        % Optimisation (only if distance_safe_flag == 1)
                        opt_time = tic;
                        if DynOpt.distance_safe_flag == 1 || DynOpt.always_opt == 1
                            
                            % save J before the optimisation to confront it later
                            if DynOpt.input_tuning == 0
                                J_before = DynOpt.cost_function(DynOpt.temp_x0,params);
                            else
                                J_before = DynOpt.cost_function(input_state,params);
                            end
                            
                            try
                                % local params
                                params_local = params;
                                
                                if DynOpt.input_tuning == 0
                                    %%%%% OPTIMISATION - NORMAL MODE %%%%%%
                                    if DynOpt.fcon_flag == 0 
                                        [NewXopt, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local),DynOpt.temp_x0,myoptioptions);
                                    else
                                        if DynOpt.multistart == 0 && DynOpt.globalsearch == 0
                                            [NewXopt, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local),DynOpt.temp_x0,...
                                                                            DynOpt.Acon, DynOpt.Bcon, DynOpt.Acon_eq, DynOpt.Bcon_eq,...
                                                                            DynOpt.lb,DynOpt.ub,DynOpt.nonlcon,myoptioptions);
                                        elseif DynOpt.multistart == 1 && DynOpt.globalsearch == 0
                                            problem = createOptimProblem('fmincon','objective',@(x)DynOpt.cost_function(x,params_local),'x0',DynOpt.temp_x0,'lb',DynOpt.lb,'ub',DynOpt.ub,'nonlcon',DynOpt.nonlcon,'options',myoptioptions);
                                            gs = MultiStart;
                                            [NewXopt, J, DynOpt.exitflag] = run(gs,problem,DynOpt.nstart); 
                                        elseif DynOpt.globalsearch == 1 && DynOpt.multistart == 0
                                            problem = createOptimProblem('fmincon','objective',@(x)DynOpt.cost_function(x,params_local),'x0',DynOpt.temp_x0,'lb',DynOpt.lb,'ub',DynOpt.ub,'nonlcon',DynOpt.nonlcon,'options',myoptioptions);
                                            gs = GlobalSearch;
                                            [NewXopt, J, DynOpt.exitflag] = run(gs,problem); 
                                        else
                                            disp('WRONG OPTIMISATION FLAGS')
                                            return
                                        end
                                    end
                                else
                                    %%% OPTIMISATION - INPUT TUNE MODE %%%%
                                    if DynOpt.fcon_flag == 0 
                                        [DynOpt.theta_u, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local),input_state,myoptioptions);
                                        NewXopt = [DynOpt.temp_x0(1:DynOpt.StateDim+DynOpt.nparams); DynOpt.theta_u; DynOpt.temp_x0(end)];
                                    else
                                        disp('WRONG OPTIMISATION FLAGS')
                                        return
                                    end
                                end
                            catch
                                params.SatellitesCoordinates = recovery_pos;
                                params.SatellitesAttitude = recovery_att;
                                NewXopt = DynOpt.temp_x0;
                                J = DynOpt.cost_function(NewXopt,params);
                            end
                            
                            % opt counter
                            DynOpt.opt_counter = DynOpt.opt_counter + 1;
                        else
                            params.SatellitesCoordinates = recovery_pos;
                            params.SatellitesAttitude = recovery_att;
                            NewXopt = DynOpt.temp_x0;
                            
                            % no care condition 
                            J_before = 1;
                            J = J_before;
                        end
                        
                        
                        % adaptive buffer backup restore
                        DynOpt.Y_space = Y_space_backup;
                        DynOpt.Y_space_full_story = Y_space_full_story_backup;
                        
                        % check J_dot condition
                        J_diff = (J/J_before);
                        distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);

                        if ( (J_diff <= DynOpt.Jdot_thresh) || (distance > DynOpt.safety_interval) )  || DynOpt.blue_flag
                            
                            % assign optimised state
                            DynOpt.X = NewXopt;
                            
                            % store measure times
                            DynOpt.opt_chosen_time = [DynOpt.opt_chosen_time k];
                            
                            % counters
                            DynOpt.jump_flag = 0;
                            DynOpt.select_counter = DynOpt.select_counter + 1;
                            
                            if (DynOpt.input_tuning == 0) || (DynOpt.optimise_input == 1)
                                DynOpt.OptXstory(:,DynOpt.BackTimeIndex) = DynOpt.X;

                                % params and state update
                                if DynOpt.identify == 1
                                    params_update(DynOpt.X);
                                end
                                x_propagate = DynOpt.X;

                                %%%%%%%%%%%%%%%%% FIRST MEASURE UPDATE %%%%%%%%
                                % manage measurements
                                % set the derivative buffer as before the optimisation process (multiple f computation)
                                back_time = DynOpt.BackTimeIndex;
                                if (back_time) >= DynOpt.d1_derivative
                                    DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
                                else
                                    init_pos = DynOpt.d1_derivative-back_time;
                                    DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
                                end

                                %%%% ESTIMATED measurements
                                % measures       
                                % NB: the output storage has to be done in
                                % back_time+1 as the propagation has been
                                % performed 
                                [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(x_propagate,0,measure_forward,DynOpt.buf_dyhat_temp,DynOpt.intYhat_full_story,params);
                                DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
                                DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
                                DynOpt.intYhat_full_story(:,back_time+1) = Yhat(:,3);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

                                %%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%
                                n_iter_propagate = DynOpt.ActualTimeIndex-DynOpt.BackTimeIndex;
                                for j =1:n_iter_propagate 
                                    % back time
                                    back_time = DynOpt.BackTimeIndex+j;

                                    % integrate
                                    [x_propagate, params] = model_propagate_local(back_time,DynOpt.Ts,x_propagate, params);                      
                                    DynOpt.OptXstory(:,back_time) = x_propagate;

                                    % manage measurements
                                    % set the derivative buffer as before the optimisation process (multiple f computation)
                                    if (back_time) >= DynOpt.d1_derivative
                                        DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
                                    else
                                        init_pos = DynOpt.d1_derivative-back_time;
                                        DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
                                    end

                                    %%%% ESTIMATED measurements
                                    % measures       
                                    % NB: the output storage has to be done in
                                    % back_time+1 as the propagation has been
                                    % performed 
                                    [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(x_propagate,0,measure_forward,DynOpt.buf_dyhat_temp,DynOpt.intYhat_full_story,params);
                                    DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
                                    DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
                                    DynOpt.intYhat_full_story(:,back_time+1) = Yhat(:,3);
                                end
                            else                               
                                % update with TRUE data
                                DynOpt.OptXstory(end-2:end-1,DynOpt.BackTimeIndex) = DynOpt.theta_u;
                                n_iter_propagate = DynOpt.ActualTimeIndex-DynOpt.BackTimeIndex;
                                for j =1:n_iter_propagate 
                                    back_time = DynOpt.BackTimeIndex+j;
                                    DynOpt.OptXstory(end-2:end-1,back_time) = DynOpt.theta_u;
                                end
                                
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                            DynOpt.Jstory(1,end+1) = J;
                        else
                            % keep the initial guess
                            DynOpt.X = DynOpt.temp_x0;
                        end                        
                        
                        % adaptive sampling
                        DynOpt.Y_space(1:end-1) = DynOpt.Y_space(2:end);
                        DynOpt.Y_space(end) = DynOpt.ActualTimeIndex;
                        DynOpt.Y_space_full_story(end+1) = DynOpt.ActualTimeIndex;

                        % stop time counter
                        DynOpt.opt_time(end+1) = toc(opt_time);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    %%% backward optimization %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else 
                        %%%% TO BE DONE %%%%
                    end
                end
                clc;
            end
        end
        DynOpt.run_time = toc(run_time);
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

