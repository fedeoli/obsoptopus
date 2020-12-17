function [DynOpt, params] = MainOpt_DEZ_general_v10_fun_params(struct)

%% Init Section
close all
clc

% dependencies
addpath(genpath([pwd '/Lib']));

% new random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% global data structure
global DynOpt params

%: Sampling Time
DynOpt.Ts = struct.Ts; 
DynOpt.Tstart = struct.T0;
DynOpt.Tend = struct.Tend;

% simulation time
DynOpt.time = DynOpt.Tstart:DynOpt.Ts:DynOpt.Tend;
DynOpt.Niter = length(DynOpt.time);

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
DynOpt.c1_dderivative = 5;
%has to be smaller than DynOpt.Nts
DynOpt.d1_dderivative = 10;

% buffer init
DynOpt.buf_dy = zeros(DynOpt.dim_out,DynOpt.d1_derivative);
DynOpt.buf_dyhat = zeros(DynOpt.dim_out,DynOpt.d1_derivative);

% gradient buffer init
DynOpt.buf_dyhat_grad = zeros(1,DynOpt.d1_derivative);

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
    state_init = DynOpt.init_error_amp.*randn(DynOpt.StateDim,1) + (1+DynOpt.init_error_amp.*randn(DynOpt.StateDim,1)).*DynOpt.Xtrue(1:DynOpt.StateDim);
    param_init = DynOpt.init_param_error_amp.*randn(length(DynOpt.param_estimate),1) + (1+DynOpt.init_param_error_amp.*randn(length(DynOpt.param_estimate),1)).*DynOpt.Xtrue(DynOpt.StateDim+1:end);
    DynOpt.X  = [state_init; param_init];
    DynOpt.X_init = DynOpt.X;
    DynOpt.Xtrue_init = DynOpt.Xtrue;
    
    % gradient descent
    DynOpt.y_end = struct.y_end;
    DynOpt.alpha_grad = struct.alpha_grad;
    DynOpt.max_iter = struct.max_iter;
    DynOpt.grad_thresh = struct.grad_thresh;
    
    % optimset
    DynOpt.get_measure = @EvaluateCostFunctionOnWindow_Output_v1;
    DynOpt.cost_function = @EvaluateCostFunctionOnWindow_general_v4;

    % filtering options
    DynOpt.filter_flag = 0;
    DynOpt.filter_window = 10;

    %%%%%%% storage %%%%%%%%
    DynOpt.WindowSamples = max(2,DynOpt.Nts*(DynOpt.w-1)+1);
    
    % measure buffer
    DynOpt.Y =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.Y_story = zeros(DynOpt.dim_out,1);
    DynOpt.Y_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.Yhat_full_story = zeros(DynOpt.dim_out,1);
    
    % measure derivative buffer
    DynOpt.dY =  zeros(DynOpt.dim_out,DynOpt.w);
    DynOpt.dY_full_story = zeros(DynOpt.dim_out,1);
    DynOpt.dYhat_full_story = zeros(DynOpt.dim_out,1);
    
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
    DynOpt.J_meas = 1;
    DynOpt.J_der = 1;
    DynOpt.temp_time = [];
    DynOpt.grad_story = zeros(DynOpt.StateDim+DynOpt.nparams,1);

    % check update condition thresholds (if blue_flag==1 the previous are not checked)
    DynOpt.J_thresh = 1e-2;
    DynOpt.diff_thresh = 1e1;
    DynOpt.blue_flag = 1;

%% Observer implementation
    disp('Processing data with the optimization-based observer...')
    tic
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
            [DynOpt.buf_dy,Y_true] = DynOpt.get_measure(DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dy);
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

        % fisrt bunch of data - read Y every Nts
        if(mod(k,DynOpt.Nts) ~= 1)
            %%%% ESTIMATED measurements
            % measures
            [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat);
            DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
            DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
            
            % clean 
            clc
        else

            % Display iteration slengthtep
            disp(['n window: ', num2str(DynOpt.w),'  n samples: ', num2str(DynOpt.Nts)])
            disp(['Iteration Number: ', num2str(floor(k/(DynOpt.Nts+1))),'/',num2str(floor(length(DynOpt.time)/(DynOpt.Nts+1)))])
            disp(['Last cost function: ', num2str(DynOpt.Jstory(end))])
            
            %%%% OUTPUT measurements - buffer of w elements
            DynOpt.Y(:,1:end-1) = DynOpt.Y(:,2:end);
            DynOpt.Y(:,end) = Y_filter(:,1);
            
            % measures derivative
            DynOpt.dY(:,1:end-1) = DynOpt.dY(:,2:end);
            DynOpt.dY(:,end) = Y_filter(:,2);
            
            % store the measure times
            DynOpt.temp_time = [DynOpt.temp_time k];

            if (k < max(1,DynOpt.WindowSamples)) 
                [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat);
                DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
            else               
                % measures
                [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat);
                DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
                DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);

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
                    
                    % store derivative buffer at this point
                    DynOpt.buf_dyhat_temp = DynOpt.buf_dyhat;

                    % Optimisation
                    NewXopt = gradient_descent_params_v2(DynOpt.xstart,DynOpt.alpha_grad,DynOpt.max_iter,DynOpt.grad_thresh);
                    J = EvaluateCostFunctionOnWindow_sensitivity(NewXopt);

                    % check J_dot condition
                    [DynOpt.J_buf, J_dot] = IterativePseudoDerivative(DynOpt.Ts*DynOpt.Nts,J,DynOpt.J_c1_derivative,DynOpt.J_d1_derivative,0,DynOpt.J_buf);
                    
                    % check state variation condition
                    state_diff = norm(DynOpt.X-NewXopt);

                    if (abs(J_dot) > DynOpt.J_thresh) || (state_diff < DynOpt.diff_thresh) || DynOpt.blue_flag
                        % store and propagate
                        DynOpt.X = NewXopt;
                        DynOpt.OptXstory(:,DynOpt.BackTimeIndex) = DynOpt.X;

                        params_update(DynOpt.X);
                        x_propagate = DynOpt.X;
                        for j =1:DynOpt.WindowSamples-1     %in absolute time as: k-(w-1)*Nts+1:k,
                            set_input(DynOpt.BackTimeIndex+j);
                            x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,1,params);
                            DynOpt.OptXstory(:,DynOpt.BackTimeIndex+j) = x_propagate;
                        end
                    else
                        NewXopt = DynOpt.temp_x0;
                        J = DynOpt.cost_function(NewXopt);
                    end
                    
                    %%%% ESTIMATED measurements
                    % update derivative buffer
                    DynOpt.buf_dyhat = DynOpt.buf_dyhat_temp;
                    % measures                      
                    [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(NewXopt,DynOpt.w,measure_forward,DynOpt.buf_dyhat_temp);
                    DynOpt.Yhat_full_story(:,end) = Yhat(:,1);
                    DynOpt.dYhat_full_story(:,end) = Yhat(:,2);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                %%% backward optimization %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else 
                    %%%% TO BE DONE %%%%
                end
                DynOpt.Jstory(1,end+1) = J;
            end
            clc;
        end
    end
    toc
    
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
    for i=1:DynOpt.y_end
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

