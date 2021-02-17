function [DynOpt, params] = MainOpt_DEZ_general_v24_fun_params(struct)

%% Init Section
close all
clc

% dependencies
addpath(genpath([pwd '/Lib']));

% set representation format
format short e

% new random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% global data structure
global DynOpt params

% identify
DynOpt.identify = struct.identify;

% raw and filtered data
if struct.filter == 1 && struct.simulationModel == 0
    DynOpt.data_raw = struct.dati_raw;
    DynOpt.data = struct.dati;
end

% Sampling Time
DynOpt.Ts = struct.Ts;
DynOpt.Tstart = struct.T0;
DynOpt.Tend = struct.Tend;

% integration method
DynOpt.int_flag = struct.int_flag;
DynOpt.tspan = [1, 1+DynOpt.Ts];

% simulation time
DynOpt.time = DynOpt.Tstart:DynOpt.Ts:DynOpt.Tend;
DynOpt.Niter = length(DynOpt.time);
DynOpt.check = struct.check;

% model name
DynOpt.modelname = struct.model;

% observer setup
DynOpt.ObserverOn = struct.ObserverOn;
DynOpt.OptimisationOn = struct.OptimisationOn;
% starts finding the state/parameters backward, -1, (actual time t_k) or Forward, 1 (i.e. at time  t_{k-(w-1)*Nts}
DynOpt.ForwardOptimization = struct.forward;

% Inter-sampling Time, number of samples within intersampled data
DynOpt.Nts = struct.Nts; 
%number of inter-sampled data in a window: the total number of sampled data in a window are (DynOpt.w-1)*DynOpt.Nts+1
DynOpt.w = struct.w;
% measure noise
DynOpt.measure_amp = struct.noise_amp;

% scale factor
DynOpt.scale_factor = struct.scale_factor;

% simulate the plant
DynOpt.simulationModel = struct.simulationModel;

% measure
DynOpt.dim_out = struct.dim_out;

% montecarlo flag
DynOpt.montecarlo = struct.montecarlo;
if DynOpt.montecarlo == 1
    % state bounds
    DynOpt.params_init = struct.params_init;
end

% init synthetic data flag
DynOpt.synthetic_flag = 0;

%% PLANT model and data
switch DynOpt.modelname
    case 'tokamak'
        tokamak_init;
        DynOpt.model = @tokamak_model_v2;
        DynOpt.param_estimate = [params.gamma0, params.gamma1];
        DynOpt.n_param_est = length(DynOpt.param_estimate);
        DynOpt.init_state = [0;0.01;0;0];
        DynOpt.StateDim = DynOpt.params.StateDim;
        
    case 'cubli'
        cubli_init;
        DynOpt.model = @cubli_model_v2_notime;
        DynOpt.param_estimate = [params.M, params.If];
        DynOpt.n_param_est = length(DynOpt.param_estimate);
        DynOpt.init_state = [0.1;0.1;0.1;0.1];
        DynOpt.StateDim = DynOpt.params.StateDim;
        
    case 'runaway'
        if DynOpt.montecarlo == 1
            runaway_init_montecarlo(DynOpt.params_init)
        else
            runaway_init_v2(struct);
        end
        DynOpt.model = @runaway_model_notime_v2;
        DynOpt.model_flow = @runaway_dflow_parameters;
        DynOpt.n_param_est = length(DynOpt.param_estimate);
        DynOpt.init_state = [params.T0; params.W0; params.x10; params.x20];
        DynOpt.StateDim = DynOpt.params.StateDim;
        
    otherwise
        disp('Model not found')
        return
end

% measures and cost function version
if DynOpt.simulationModel == 1
    DynOpt.get_measure = @get_measure_v5;
    DynOpt.get_measure_name = 'get_measure_v5';
else
    DynOpt.get_measure = @get_measure_v5;
    DynOpt.get_measure_name = 'get_measure_v5'; 
    DynOpt.get_measure_data = @get_measure_v5_data;
    DynOpt.get_measure_name_data = 'get_measure_v5_data'; 
end
DynOpt.cost_function = @cost_function_v4;
DynOpt.cost_function_name = 'cost_function_v4';

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
DynOpt.J_c1_derivative = 3;
DynOpt.J_d1_derivative = 10;
DynOpt.J_buf = zeros(1,DynOpt.J_d1_derivative);
    
%% model simulation
tic
if(DynOpt.simulationModel == 1)
    
    % display sim data
    disp('Model simulation')
    disp(['gamma: ', num2str(params.gamma), ' gamma1: ', num2str(params.gamma1), ' Wt: ', num2str(params.Wt)])
    
    % notify that you're generating synthetic data
    DynOpt.synthetic_flag = 1;
    
    % state storage
    DynOpt.stateStory = zeros(DynOpt.StateDim,DynOpt.Niter);
    DynOpt.stateStory(:,1) = DynOpt.init_state;
    
    % temp vars
    junk_buf = zeros(DynOpt.dim_out,DynOpt.d1_derivative);
    
    % output storage
    DynOpt.outputStory = zeros(DynOpt.dim_out,DynOpt.Niter);
    [~,temp_y] = DynOpt.get_measure(DynOpt.init_state,0,1,junk_buf,junk_buf,params);
    DynOpt.outputStory(1) = temp_y(1);
         
    % model integration
    for k=2:DynOpt.Niter  
        
        % current time
        DynOpt.ActualTimeIndex = k;
        
        % system control input
        set_input(k);
        
        % integration
        temp_state = rk4_V1_1(DynOpt.model, DynOpt.tspan, DynOpt.stateStory(:,k-1), params);   
        
        % state storage
        DynOpt.stateStory(:,k) = temp_state(:,end);
        [~, temp_y] = DynOpt.get_measure(DynOpt.stateStory(:,k),0,1,junk_buf,junk_buf,params);
        DynOpt.outputStory(k) = temp_y(1);
    end     
    
    % stop generating synthetic data
    DynOpt.synthetic_flag = 0;
end
toc


%% OBSERVER SETUP

if DynOpt.ObserverOn == 1
    %true state and parameters
    if DynOpt.simulationModel == 1
        DynOpt.Xtrue = [DynOpt.stateStory(:,1); DynOpt.param_estimate'];
    else
        DynOpt.Xtrue = [DynOpt.params.T0; DynOpt.params.W0; DynOpt.params.x10; DynOpt.params.x20; DynOpt.param_estimate'];
    end
    
    % state length
    DynOpt.aug_state_dim = length(DynOpt.Xtrue);
    DynOpt.nparams = length(DynOpt.param_estimate);
    
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
    DynOpt.else_flag = 0;
    DynOpt.safety_density = struct.safety_density;
    DynOpt.safety_interval = int32(DynOpt.safety_density*DynOpt.WindowSamples);
    
    % dJ condition buffer (adaptive sampling)
    DynOpt.dJ_cond_story = zeros(3,1);
    
    % performance analysis
    DynOpt.performance_n_rows = 2*DynOpt.Nts;
    DynOpt.lambda_min = 0;
    DynOpt.measure_exp = 1;

    % corrupted trajectory
    DynOpt.Xstory = zeros(length(DynOpt.X),length(DynOpt.time));
    DynOpt.Xstory(:,1) = DynOpt.X_init;

    % optimisation trajectory
    DynOpt.OptXstory = zeros(length(DynOpt.X),length(DynOpt.time));
    DynOpt.OptXstory(:,1) = DynOpt.X_init;

    % optimisation error
    DynOpt.OptErrorStory = DynOpt.OptXstory;
    
    % set reference trajectory
    if DynOpt.simulationModel == 1
        DynOpt.OptXstoryTRUE = [DynOpt.stateStory(1:DynOpt.StateDim,:);DynOpt.param_estimate'.*ones(length(DynOpt.param_estimate),length(DynOpt.time))];
    else
        DynOpt.OptXstoryTRUE = DynOpt.dataStory;
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
    
    % safety interval policy
    DynOpt.always_opt = struct.always_opt;
    
    % opt counters
    DynOpt.opt_counter = 0;
    DynOpt.select_counter = 0;

%% Optimisation implementation
    if DynOpt.OptimisationOn
        disp('Processing data with the optimization-based observer...')
        run_time = tic;
        for k=1:length(DynOpt.time)

            % update actual index
            DynOpt.ActualTimeIndex = k;
            
            % true state
            if DynOpt.simulationModel == 1
                DynOpt.Xtrue = [DynOpt.stateStory(:,k);DynOpt.param_estimate'];
            else
                DynOpt.Xtrue = DynOpt.dataStory(k);
            end

            %forward propagation of the previous estimate
            if(k>1)

                % input
                set_input(k);

                % integration 
                temp = rk4_V1_1(DynOpt.model, DynOpt.tspan, DynOpt.OptXstory(:,k-1), params);
                DynOpt.X = temp(:,end);
                DynOpt.OptXstory(:,k) = DynOpt.X; 
                
                % integration - without opt
                temp = rk4_V1_1(DynOpt.model, DynOpt.tspan, DynOpt.Xstory(:,k-1), params);
                DynOpt.Xstory(:,k) = temp(:,end);
                           
            end

            %%%%%%%%%% MEASUREMENT %%%%%%%%%%%
            % read measure 
            measure_forward = 1;

            if DynOpt.simulationModel == 1
                [DynOpt.buf_dy,Y_true] = DynOpt.get_measure(DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dy,DynOpt.intY_full_story,params);
                DynOpt.measure_noise = DynOpt.measure_amp*randn(1);
                % copy to Y and add noise 
                Y_noise = Y_true + [DynOpt.measure_noise; zeros(DynOpt.dim_out-1,1)];
            else
                [DynOpt.buf_dy,Y_true] = DynOpt.get_measure_data(DynOpt.ActualTimeIndex,measure_forward,DynOpt.buf_dy,DynOpt.dataStory);
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
            % rich enough by using the derivative
            dJ_cond(DynOpt.theta);
            distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);
            DynOpt.distance_safe_flag = (distance < DynOpt.safety_interval);
            % SENSOR SATURATION/ACCURACY
            DynOpt.sensor_stop_flag = (DynOpt.Y_full_story(DynOpt.ActualTimeIndex) < DynOpt.sensor_stop_thresh);
            % if last sample is too near or the condition isn't met then
            % propagate. However, check that the last optimisation has been
            % performed at most WindowSamples times ago
            if  (((distance < DynOpt.Nts) || (DynOpt.dJ_cond < DynOpt.dJcond_thresh)) && DynOpt.distance_safe_flag ) || (DynOpt.sensor_stop_flag)
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
                if DynOpt.simulationModel == 1
                    disp(['TARGET: gamma: ', num2str(DynOpt.params.gamma), ' gamma1: ', num2str(DynOpt.params.gamma1)])
                end
                disp(['INIT: gamma: ', num2str(DynOpt.X_init(end-1)), ' gamma1: ', num2str(DynOpt.X_init(end))])
                disp(['CURRENT: gamma: ', num2str(params.gamma), ' gamma1: ', num2str(params.gamma1)])
                disp(['n window: ', num2str(DynOpt.w),'  n samples: ', num2str(DynOpt.Nts)])
                disp(['Iteration Number: ', num2str(k),'/',num2str(length(DynOpt.time))])
                disp(['Last cost function: ', num2str(DynOpt.Jstory(end))]);
                disp(['N. optimisations RUN: ',num2str(DynOpt.opt_counter)]);
                disp(['N. optimisations SELECTED: ',num2str(DynOpt.select_counter)]);

                %%%% OUTPUT measurements - buffer of w elements
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

                % store the measure times
                DynOpt.temp_time = [DynOpt.temp_time k];

                % wait until enough samples have been collected
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

                        % Optimisation
                        opt_time = tic;
                        if DynOpt.distance_safe_flag == 1 || DynOpt.always_opt == 1
                            
                            % save J before the optimisation to confront it later
                            J_before = DynOpt.cost_function(DynOpt.temp_x0,params);
                            
                            try
                                % local params
                                params_local = params;

                                % optimisation
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
                            catch
                                NewXopt = DynOpt.temp_x0;
                                J = DynOpt.cost_function(NewXopt,params);
                            end

                            % opt counter
                            DynOpt.opt_counter = DynOpt.opt_counter + 1;
                        else
                            % no care condition 
                            NewXopt = DynOpt.temp_x0;
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
                            % setup
                            DynOpt.X = NewXopt;
                            DynOpt.OptXstory(:,DynOpt.BackTimeIndex) = DynOpt.X;
                            
                            % store measure times
                            DynOpt.opt_chosen_time = [DynOpt.opt_chosen_time k];

                            % counters
                            DynOpt.jump_flag = 0;
                            DynOpt.select_counter = DynOpt.select_counter + 1;

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
                            DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                            %%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%
                            n_iter_propagate = DynOpt.ActualTimeIndex-DynOpt.BackTimeIndex;
                            for j=1:n_iter_propagate     
                                % back time
                                back_time = DynOpt.BackTimeIndex+j;

                                % set input
                                set_input(back_time);

                                % integrate
                                temp = rk4_V1_1(DynOpt.model, DynOpt.tspan, x_propagate, params);
                                x_propagate = temp(:,end);                         
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
                                DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);
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
        DynOpt.Ytrue_full_story = DynOpt.Ytrue_full_story(:,2:end);
        DynOpt.dY_full_story = DynOpt.dY_full_story(:,2:end);
        DynOpt.dYhat_full_story = DynOpt.dYhat_full_story(:,2:end);
        DynOpt.Jstory = DynOpt.Jstory(:,2:end);
        DynOpt.Y_space_full_story = diff(DynOpt.Y_space_full_story);
        DynOpt.intY_full_story = DynOpt.intY_full_story(:,2:end);
        DynOpt.intYhat_full_story = DynOpt.intYhat_full_story(:,2:end);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    else
        % display sim data
        disp('Model simulation')
        disp(['gamma: ', num2str(params.gamma), ' gamma1: ', num2str(params.gamma1), ' Wt: ', num2str(params.Wt)])

        % model integration - ask for stop flag settings/meaning
        for k=2:DynOpt.Niter       
            % system control input
            set_input(k);
            temp_state = rk4_V1_1(DynOpt.model, DynOpt.tspan, DynOpt.OptXstory(:,k-1), params);   
            DynOpt.OptXstory(:,k) = temp_state(:,end);
        end
    end
else
    %% TBD
end

end

