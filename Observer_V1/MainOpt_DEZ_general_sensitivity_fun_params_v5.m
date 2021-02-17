function [DynOpt, params] = MainOpt_DEZ_general_sensitivity_fun_params_v5(struct)

%% Init Section
close all
clc

% dependencies
addpath(genpath([pwd '/Lib']));

% new random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% global data structure
global DynOpt params

% raw and filtered data
if struct.filter == 1 && struct.simulationModel == 0
    DynOpt.data_raw = struct.dati_raw;
    DynOpt.data = struct.dati;
end

%: Sampling Time
DynOpt.Ts = struct.Ts; 
DynOpt.Tstart = struct.T0;
DynOpt.Tend = struct.Tend;

% simulation time
DynOpt.time = DynOpt.Tstart:DynOpt.Ts:DynOpt.Tend;
DynOpt.Niter = length(DynOpt.time);
DynOpt.check = struct.check;

% model name
DynOpt.modelname = struct.model;

% observer setup
DynOpt.ObserverOn = struct.ObserverOn;
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
DynOpt.check = struct.check;

% measure
DynOpt.dim_out = struct.dim_out;


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
        DynOpt.model = @cubli_model_v2;
        DynOpt.param_estimate = [params.M, params.If];
        DynOpt.n_param_est = length(DynOpt.param_estimate);
        DynOpt.init_state = [0.1;0.1;0.1;0.1];
        DynOpt.StateDim = DynOpt.params.StateDim;
        
    case 'runaway'
        runaway_init;
        DynOpt.model = @runaway_model;
        DynOpt.model_flow = @runaway_dflow_parameters;
        DynOpt.param_estimate = [params.gamma, params.gamma1];
        DynOpt.n_param_est = length(DynOpt.param_estimate);
        DynOpt.init_state = [params.T0; params.W0];
        DynOpt.StateDim = DynOpt.params.StateDim;
        
    otherwise
        disp('Model not found')
        return
end
    
% model simulation
tic
if(DynOpt.simulationModel == 1)
    
    disp('Model simulation')
    
    % state storage
    DynOpt.stateStory = zeros(DynOpt.StateDim,DynOpt.Niter);
    DynOpt.stateStory(:,1) = DynOpt.init_state;
    
    % model integration - ask for stop flag settings/meaning
    for k=2:DynOpt.Niter       
        % system control input
        set_input(k);
        DynOpt.stateStory(:,k) = PlantJumpMap_general_notime_params(DynOpt.stateStory(:,k-1),DynOpt.model,1,params);            
    end     
else
    disp('Get data')
    DynOpt.Tstart = struct.T0;
    DynOpt.Tend = struct.Tend;
    DynOpt.sample_time = struct.sample_time;
    
    [rows,cols] = size(struct.dati);
    
    [junk,pos_0] = min(abs(struct.dati(:,1)-struct.T0));
    [junk,pos_end] = min(abs(struct.dati(:,1)-struct.Tend));
    
    time_sample = pos_0:DynOpt.sample_time:pos_end;
    values = struct.dati(time_sample,2);
    
    DynOpt.stateStory = zeros(cols-1,length(values));
    DynOpt.stateStory(1,:) = values;
    DynOpt.time = struct.dati(time_sample);
end
toc
%% Derivatives setup

disp('Setting derivatives')

% first derivative 
% number of derivatives
DynOpt.c1_derivative = 1;
%has to be smaller than DynOpt.Nts
DynOpt.d1_derivative = 3;

% second derivative
DynOpt.c1_dderivative = 1;
%has to be smaller than DynOpt.Nts
DynOpt.d1_dderivative = 3;

% buffer init
DynOpt.buf_dy = zeros(DynOpt.dim_out,DynOpt.d1_derivative);
DynOpt.buf_dyhat = zeros(DynOpt.dim_out,DynOpt.d1_derivative);
DynOpt.buf_ddy = zeros(DynOpt.dim_out,DynOpt.d1_dderivative);
DynOpt.buf_ddyhat = zeros(DynOpt.dim_out,DynOpt.d1_dderivative);

% gradient buffer init
% DynOpt.buf_dyhat_grad = zeros(1,DynOpt.d1_derivative);
% DynOpt.buf_ddyhat_grad = zeros(1,DynOpt.d1_dderivative);

% cost function derivative
DynOpt.J_c1_derivative = 6;
DynOpt.J_d1_derivative = 20;
DynOpt.J_buf = zeros(1,DynOpt.J_d1_derivative);

if(DynOpt.simulationModel == 1)
    DynOpt.state = DynOpt.stateStory;
else
    DynOpt.state = DynOpt.stateStory;
end


%% OBSERVER SETUP

if DynOpt.ObserverOn == 1

    % nparams
    DynOpt.nparams = length(DynOpt.param_estimate);

    %true state and parameters
    if DynOpt.simulationModel == 1
        DynOpt.Xtrue = [DynOpt.state(:,1); DynOpt.param_estimate'];
    else
        DynOpt.Xtrue = [DynOpt.params.T0; DynOpt.params.W0; DynOpt.param_estimate'];
    end

    DynOpt.aug_state_dim = length(DynOpt.Xtrue);
    DynOpt.init_error_amp = struct.init_error_amp;
    DynOpt.init_param_error_amp = struct.init_param_error_amp;
    DynOpt.paramsAllPos = struct.ParamsAllPos;
    state_init = DynOpt.init_error_amp.*randn(DynOpt.StateDim,1) + (1+DynOpt.init_error_amp.*randn(DynOpt.StateDim,1)).*DynOpt.Xtrue(1:DynOpt.StateDim);
    param_init = DynOpt.init_param_error_amp.*randn(length(DynOpt.param_estimate),1) + (1+DynOpt.init_param_error_amp.*randn(length(DynOpt.param_estimate),1)).*DynOpt.Xtrue(DynOpt.StateDim+1:end);
    if DynOpt.paramsAllPos
        param_init = abs(param_init);
    end
    DynOpt.X  = [state_init; param_init];
    DynOpt.X_init = DynOpt.X;
    DynOpt.Xtrue_init = DynOpt.Xtrue;
    
    % gradient descent
    DynOpt.y_end = struct.y_end;
    DynOpt.alpha_grad = struct.alpha_grad;
    DynOpt.max_iter = struct.max_iter;
    DynOpt.identify = struct.identify;
    DynOpt.grad_thresh = struct.grad_thresh;
    DynOpt.opt_time = 0;
    DynOpt.alpha_dyn = struct.alpha_dyn;
    
    % optimset
    DynOpt.get_measure = @EvaluateCostFunctionOnWindow_Output_sensitivity_v1;
    DynOpt.cost_function = @EvaluateCostFunctionOnWindow_general_v3;

    % filtering options
    DynOpt.filter_flag = 0;
    DynOpt.filter_window = 10;

    %%%%%%% storage %%%%%%%%
    DynOpt.WindowSamples = max(2,DynOpt.Nts*(DynOpt.w-1)+1);
    
    % measure buffer
    DynOpt.Y =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.Y_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.Yhat_full_story = zeros(DynOpt.dim_out,1);
    
    % measure derivative buffer
    DynOpt.dY =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.dY_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.dYhat_full_story = zeros(DynOpt.dim_out,1);
    
    % measure second derivative buffer
    DynOpt.ddY =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.ddY_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.ddYhat_full_story = zeros(DynOpt.dim_out,1);
    
    % buffer adaptive sampling
    DynOpt.Y_space = zeros(1,DynOpt.w);
    DynOpt.Y_space_full_story = 0;
    DynOpt.dJcond_thresh = struct.dJcond_thresh;
    DynOpt.theta = struct.theta;
    DynOpt.else_flag = 0;
    
    % dJ condition buffer (adaptive sampling)
    DynOpt.dJ_cond_story = zeros(3,1);
    
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
        DynOpt.OptXstoryTRUE = [DynOpt.state(1:DynOpt.StateDim,:);DynOpt.param_estimate'.*ones(2,length(DynOpt.time))];
    else
        DynOpt.OptXstoryTRUE = [DynOpt.state(1,:);DynOpt.param_estimate'.*ones(2,length(DynOpt.time))];
    end
    
    % observer cost function init
    DynOpt.J = 1E3;
    DynOpt.Jstory = DynOpt.J;
    DynOpt.Jdot_story = 0;
    DynOpt.J_meas = 1;
    DynOpt.J_der = 1;
    DynOpt.temp_time = [];
    DynOpt.grad_story = zeros(DynOpt.StateDim+DynOpt.nparams,1);

    % check update condition thresholds (if blue_flag==1 the previous are not checked)
    DynOpt.J_thresh = struct.J_thresh;
    DynOpt.Jdot_thresh = struct.Jdot_thresh;
    DynOpt.diff_thresh = 0;
    DynOpt.blue_flag = struct.blue_flag;
    
    % opt counters
    DynOpt.opt_counter = 0;
    DynOpt.select_counter = 0;

%% Observer implementation
    disp('Processing data with the optimization-based observer...')
    run_time = tic;
    for k=1:length(DynOpt.time)
        
        % update actual index
        DynOpt.ActualTimeIndex = k;
        
        %updating the moving window data with the inter-sampled data
        DynOpt.Xtrue = [DynOpt.state(:,k);DynOpt.param_estimate'];

        %forward propagation of the previous estimate
        if(k>1)
            % update parameters
            params_update(DynOpt.X);
            
            % input
            set_input(k);

            % integration
            DynOpt.X = PlantJumpMap_general_notime_params(DynOpt.OptXstory(:,k-1), DynOpt.model, 1 ,params);
            DynOpt.OptXstory(:,k) = DynOpt.X;
            DynOpt.Xstory(:,k) = PlantJumpMap_general_notime_params(DynOpt.Xstory(:,k-1), DynOpt.model ,1 , DynOpt.params);
        end
        
        %%%%%%%%%% MEASUREMENT %%%%%%%%%%%
        % read measure 
        measure_forward = 1;

        if DynOpt.simulationModel == 1
            [DynOpt.buf_dy,DynOpt.buf_ddy,Y_true] = DynOpt.get_measure(DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dy,DynOpt.buf_ddy);
            DynOpt.measure_noise = DynOpt.measure_amp*randn(1);
            % copy to Y and add noise 
            Y_noise = Y_true + [DynOpt.measure_noise; zeros(DynOpt.dim_out-1,1)];
        else
            [DynOpt.buf_dy,DynOpt.buf_ddy,DynOpt.buf_yint, Y_true] = EvaluateCostFunctionOnWindow_Output_general_data(k,measure_forward,DynOpt.buf_dy,DynOpt.buf_ddy,DynOpt.buf_yint);
            Y_noise = Y_true;
        end

        % no filtering
        Y_filter = Y_noise;
        
        % store total memory
        DynOpt.Y_full_story(:,end+1) = Y_filter(:,1);
        DynOpt.dY_full_story(:,end+1) = Y_filter(:,2);
        DynOpt.ddY_full_story(:,end+1) = Y_filter(:,3);

        % fisrt bunch of data - read Y every Nts and check if the signal is
        % rich enough by using the derivative
        dJ_cond_sensitivity_v1(DynOpt.theta);
        distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);
        if  (distance < DynOpt.Nts) || (DynOpt.dJ_cond)
            %%%% ESTIMATED measurements
            % measures
            [DynOpt.buf_dyhat, DynOpt.buf_ddyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.buf_ddyhat);
            DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
            DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
            DynOpt.ddYhat_full_story(:,end+1) = Yhat(:,3);
            
            % clean 
            clc
        else

            % Display iteration slengthtep
            disp(['TARGET: gamma: ', num2str(DynOpt.params.gamma), ' gamma1: ', num2str(DynOpt.params.gamma1)])
            disp(['INIT: gamma: ', num2str(DynOpt.X_init(3)), ' gamma1: ', num2str(DynOpt.X_init(4))])
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
            DynOpt.ddY(:,1:end-1) = DynOpt.ddY(:,2:end);
            DynOpt.ddY(:,end) = Y_filter(:,3);
            
            % store the measure times
            DynOpt.temp_time = [DynOpt.temp_time k];

            % wait until enough samples have been collected
            if (k < max(1,DynOpt.WindowSamples)) 
                [DynOpt.buf_dyhat, DynOpt.buf_ddyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.buf_ddyhat);
                DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
                DynOpt.ddYhat_full_story(:,end+1) = Yhat(:,3);
            else       
                
                % measures
                [DynOpt.buf_dyhat, DynOpt.buf_ddyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.buf_ddyhat);
                DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
                DynOpt.ddYhat_full_story(:,end+1) = Yhat(:,3);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                %%% forward optimization %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(DynOpt.ForwardOptimization == 1) 

                    % start from here - optimisation
                    DynOpt.BackTimeIndex = k-(DynOpt.w-1)*DynOpt.Nts; 

                    % set of initial conditions
                    DynOpt.temp_x0 = DynOpt.OptXstory(:,DynOpt.BackTimeIndex);
                    DynOpt.xstart.b0 = DynOpt.temp_x0;
                    DynOpt.xstart.e0 = eye(DynOpt.StateDim+DynOpt.nparams);

                    % Optimisation
                    opt_time = tic;
                    NewXopt = gradient_descent_params_v4(DynOpt.xstart,DynOpt.max_iter,DynOpt.grad_thresh);
                    J = DynOpt.cost_function(NewXopt);
                    
                    % opt counter
                    DynOpt.opt_counter = DynOpt.opt_counter + 1;

                    % check J_dot condition
                    last_Nts = DynOpt.ActualTimeIndex - DynOpt.Y_space(end);
                    J_dot = (J -DynOpt.Jstory(end))/last_Nts;

                    if (J_dot < DynOpt.Jdot_thresh) || DynOpt.blue_flag
                        % store and propagate
                        DynOpt.X = NewXopt;
                        DynOpt.OptXstory(:,DynOpt.BackTimeIndex) = DynOpt.X;
                        
                        % counters
                        DynOpt.jump_flag = 0;
                        DynOpt.select_counter = DynOpt.select_counter + 1;
                        
                        if ((DynOpt.ActualTimeIndex-DynOpt.Y_space(end)) < DynOpt.WindowSamples) || (isempty(find(DynOpt.Y_space,1))) || (DynOpt.else_flag == 1)
%                         if 1
                            DynOpt.else_flag = 0;
                            DynOpt.Y_space(1:end-1) = DynOpt.Y_space(2:end);
                            DynOpt.Y_space(end) = DynOpt.ActualTimeIndex;
                        else
                            DynOpt.else_flag = 1;
                            buf = 1;
                            DynOpt.Y_space = [DynOpt.Y_space(end-buf+1):DynOpt.Y_space(end), zeros(1,DynOpt.w-buf)];
                            DynOpt.Y = [DynOpt.Y(end-buf+1):DynOpt.Y(end), zeros(1,DynOpt.w-buf)];
                            DynOpt.dY = [DynOpt.dY(end-buf+1):DynOpt.dY(end), zeros(1,DynOpt.w-buf)];
                            DynOpt.ddY = [DynOpt.ddY(end-buf+1):DynOpt.ddY(end), zeros(1,DynOpt.w-buf)];
                        end
                        DynOpt.Y_space_full_story(end+1) = DynOpt.ActualTimeIndex;

                        % params and state update
                        if DynOpt.identify == 1
                            params_update(DynOpt.X);
                        end
                        x_propagate = DynOpt.X;
                        
                        % propagate state (x_propagate is used inside the script)
                        state_update;
                        
                        % compute and store J and Jdot
                        last_Nts = DynOpt.Y_space(end) - DynOpt.Y_space(end-1);
                        J_dot = (J -DynOpt.Jstory(end))/last_Nts;
                        DynOpt.Jdot_story(1,end+1) = J_dot;
                        DynOpt.Jstory(1,end+1) = J;
                    else
                        % keep the initial guess
                        DynOpt.X = DynOpt.temp_x0;
                        % restore the initial guess parameters, because
                        % they have been updated during the gradient
                        % descent.
                        params_update(DynOpt.X);
                    end
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
    DynOpt.dY_full_story = DynOpt.dY_full_story(:,2:end);
    DynOpt.dYhat_full_story = DynOpt.dYhat_full_story(:,2:end);
    DynOpt.Jstory = DynOpt.Jstory(:,2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % integrate derivatives to see if they're correct (only for the gyro-3 dim) 
    DynOpt.Y_int = zeros(DynOpt.y_end,DynOpt.Niter);
    DynOpt.Yhat_int = zeros(DynOpt.y_end,DynOpt.Niter);
    for i=1:1
       DynOpt.Y_int(i,:) = cumtrapz(DynOpt.Ts,DynOpt.dY_full_story(i,:)); 
       DynOpt.Yhat_int(i,:) = cumtrapz(DynOpt.Ts,DynOpt.dYhat_full_story(i,:));
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sensitivity 
    temp_time = DynOpt.temp_time;
    temp_time(find(DynOpt.temp_time<DynOpt.WindowSamples))=[];
    DynOpt.sens = zeros(length(DynOpt.X_init),length(temp_time));
    for i=1:length(DynOpt.X_init)
        DynOpt.sens(i,:) = DynOpt.Jstory./DynOpt.OptXstory(i,temp_time);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    % performance
    DynOpt.lambda_min = DynOpt.lambda_min(2:end);
    DynOpt.Jstory = DynOpt.Jstory(1,2:end);
    DynOpt.grad_story = DynOpt.grad_story(:,2:end);
    
    % Optimisation time
    DynOpt.opt_time = DynOpt.opt_time(2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

