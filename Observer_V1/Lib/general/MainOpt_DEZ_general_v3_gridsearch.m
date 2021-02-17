
%% Init Section
close all
keep counter_sample counter_window final_path
clc

% dependencies
addpath(genpath([pwd '/Lib']));

% global data structure
global DynOpt params

%: Sampling Time
DynOpt.Ts = 5e-4; 
DynOpt.Tstart = 0;
DynOpt.Tend = 1;

% simulation time
DynOpt.time = DynOpt.Tstart:DynOpt.Ts:DynOpt.Tend;
DynOpt.Niter = length(DynOpt.time);

% observer setup
% different from zero to test only the backward integration and see if it is correct: please set the correct value in the optimized variable initialization
DynOpt.TestDynamics = 0; 
% starts finding the state/parameters backward, -1, (actual time t_k) or Forward, 1 (i.e. at time  t_{k-(w-1)*Nts}
DynOpt.ForwardOptimization = 1; 

% Inter-sampling Time, number of samples within intersampled data
DynOpt.Nts = counter_sample; 
%number of inter-sampled data in a window: the total number of sampled data in a window are (DynOpt.w-1)*DynOpt.Nts+1
DynOpt.w = counter_window;
% measure noise
DynOpt.measure_amp = 5e-5;

% window weight size, last on the more recent t_k
DynOpt.Weight  = 1*ones(1,DynOpt.w); 

% window weight size, last on the more recent t_k
DynOpt.dWeight  = 0*ones(1,DynOpt.w); 
DynOpt.dWeight(1) = 0; 
DynOpt.dWeight(end) = 0;

%window weight size, last on the more recent t_k
DynOpt.ddWeight  = 0*ones(1,DynOpt.w); 
DynOpt.ddWeight(1) = 0; 
DynOpt.ddWeight(end) = 0;

% simulate the plant
simulationModel = 1;


%% PLANT model and data

model = 'tokamak';
tokamak_init;
DynOpt.model = @tokamak_model;
DynOpt.param_estimate = [params.gamma0, params.gamma1];
DynOpt.init_state = [0;0.01;0;0];

% model = 'cubli';
% cubli_init;
% DynOpt.model = @cubli_model;
% DynOpt.param_estimate = [params.M, params.If];
% DynOpt.init_state = [0.1;0.1;0.1;0.1];

    
% model simulation
if(simulationModel == 1)
    
    disp('Model simulation')
    
    % state storage
    DynOpt.stateStory = zeros(4,DynOpt.Niter);
    DynOpt.stateStory(:,1) = DynOpt.init_state;
    
    
    % model integration - ask for stop flag settings/meaning
    for k=2:DynOpt.Niter       
        % system control input
        if DynOpt.params.input_flag == 1
              params.u = params.U(k);
%             DynOpt.u = DynOpt.Km*sin(DynOpt.stateStory(3,k-1));
        else
              params.u = 0;
        end
        DynOpt.stateStory(:,k) = PlantJumpMap_general_notime(DynOpt.stateStory(:,k-1),DynOpt.model,1);            
    end          


end

%% Derivatives setup

disp('Setting derivatives')

% state selection
state_select = 1;

% first derivative 
% number of derivatives
DynOpt.c1_derivative = 6;
%has to be smaller than DynOpt.Nts
DynOpt.d1_derivative = 20;

% second derivative
DynOpt.c1_dderivative = 8;
%has to be smaller than DynOpt.Nts
DynOpt.d1_dderivative = 30;

% buffer init
DynOpt.buf_dy = zeros(1,DynOpt.d1_derivative);
DynOpt.buf_dyhat = zeros(1,DynOpt.d1_derivative);
DynOpt.buf_ddy = zeros(1,DynOpt.d1_dderivative);
DynOpt.buf_ddyhat = zeros(1,DynOpt.d1_dderivative);

if(simulationModel == 1)
    DynOpt.state = DynOpt.stateStory;
else
%%%%%%%%%%%%%%% MISSING %%%%%%%%%%%%%%%
end


%% OBSERVER SETUP

%true state and parameters
DynOpt.Xtrue = [DynOpt.state(:,1); DynOpt.param_estimate']; 

DynOpt.aug_state_dim = length(DynOpt.Xtrue);
DynOpt.init_error_amp = 2e-1*abs(mean(DynOpt.Xtrue));
DynOpt.X  = DynOpt.init_error_amp*randn(DynOpt.aug_state_dim,1) + DynOpt.Xtrue;
% DynOpt.X  = 1.2*DynOpt.Xtrue;

DynOpt.X_init = DynOpt.X;
DynOpt.Xtrue_init = DynOpt.Xtrue;

%%%%%%% storage %%%%%%%%
DynOpt.WindowSamples = max(2,DynOpt.Nts*(DynOpt.w-1)+1);
DynOpt.dim_out = 3;
DynOpt.Y =  zeros(DynOpt.dim_out,DynOpt.w);
DynOpt.Y_story = zeros(DynOpt.dim_out,1);
DynOpt.Y_full_story = zeros(DynOpt.dim_out,1);

DynOpt.performance_n_rows = 2*DynOpt.Nts;
DynOpt.lambda_min = 0;

DynOpt.Xstory = zeros(length(DynOpt.X),length(DynOpt.time));
DynOpt.Xstory(:,1) = DynOpt.X;

DynOpt.OptXstory = zeros(length(DynOpt.X),length(DynOpt.time));
DynOpt.OptXstory(:,1) = DynOpt.X;

DynOpt.OptErrorStory = DynOpt.OptXstory;
DynOpt.OptXstoryTRUE = [DynOpt.state(1:4,:);DynOpt.param_estimate'.*ones(2,length(DynOpt.time))];
DynOpt.Jstory = zeros(1,length(DynOpt.time));

% observer optimisation options
myoptioptions = optimoptions(@fminunc,'Algorithm','quasi-newton',  'MaxIter', 50,'display','off'); % 'Algorithm', 'trust-region','SpecifyObjectiveGradient', true,
%myoptioptions = optimoptions(@fmincon,  'MaxIter', 100,'display','off'); % 'Algorithm', 'trust-region','SpecifyObjectiveGradient', true,

% observer cost function init
DynOpt.J = 1E3;
temp_time = [];

if(DynOpt.TestDynamics ~= 0 )
%%%%%%%%%%%%%%% MISSING %%%%%%%%%%%%%%%
else
    
    disp('Processing data with the optimization-based observer...')
    tic
    for k=1:length(DynOpt.time)
        
        %updating the moving window data with the inter-sampled data
        DynOpt.ActualTimeIndex = k;
        DynOpt.Xtrue = [DynOpt.state(:,k);DynOpt.param_estimate'];
        
        % system control input
        if DynOpt.params.input_flag == 1
            params.u = params.U(k);
%             DynOpt.u = DynOpt.Km*sin(DynOpt.X(3));
        else
            params.u = 0;
        end
        
        if(k>1)%forward propagation of the previous estimate
            DynOpt.X = PlantJumpMap_general_notime(DynOpt.X, @cubli_model ,1);
            DynOpt.OptXstory(:,k) = DynOpt.X;
            
            DynOpt.Xstory(:,k) = PlantJumpMap_general_notime(DynOpt.Xstory(:,k-1), @cubli_model ,1);
        end
        
        % read measure 
        [DynOpt.buf_dy,DynOpt.buf_ddy,Y_true] = EvaluateCostFunctionOnWindow_Output_general_notime_nobuf(DynOpt.state(:,k),0,1,DynOpt.buf_dy,DynOpt.buf_ddy);
        DynOpt.measure_noise = DynOpt.measure_amp*randn(1);
        Y_noise = Y_true + [DynOpt.measure_noise; zeros(DynOpt.dim_out-1,1)];
        DynOpt.Y_full_story(:,end+1) = Y_noise;
        
        % fisrt bunch of data - read Y every Nts
        if(mod(k,DynOpt.Nts) == 0)
            
            % Display iteration step
            disp(['n window: ', num2str(counter_window),'  n samples: ', num2str(counter_sample)])
            disp(['Iteration Number: ', num2str(floor(k/(DynOpt.Nts+1))),'/',num2str(floor(length(DynOpt.time)/(DynOpt.Nts+1)))])
            
            %%%% OUTPUT measurements - buffer of w elements
            DynOpt.Y(:,1:end-1) = DynOpt.Y(:,2:end);
            DynOpt.Y(:,end) = Y_noise;
            
            DynOpt.Y_story(:,end+1) = Y_noise;
            temp_time = [temp_time k];
            
            % performance index computation
            if k >= (DynOpt.Nts+DynOpt.performance_n_rows)
                [DynOpt.lambda_min(end+1),~] = performance_index(DynOpt.Y_full_story,DynOpt.performance_n_rows);
            end
            
            if( k >= max(1,DynOpt.WindowSamples) )%enough samples have been acquired
                
                
                % forward optimization
                if(DynOpt.ForwardOptimization == 1) 
                    
                    DynOpt.BackTimeIndex = k-(DynOpt.w-1)*DynOpt.Nts; 
                    DynOpt.temp_x0 = DynOpt.OptXstory(:,DynOpt.BackTimeIndex);
                    [NewXopt, J] = fminunc(@EvaluateCostFunctionOnWindow_general_notime_nobuf,DynOpt.temp_x0,myoptioptions);
                    
                    DynOpt.X = NewXopt;
                    DynOpt.OptXstory(:,DynOpt.BackTimeIndex) = DynOpt.X;
                    
                    % params update
                    if strcmp(model,'cubli')
                        params.M = DynOpt.X(5);
                        params.If = DynOpt.X(6);
                    elseif strcmp(model,'tokamak')
                        params.gamma0 = DynOpt.X(5);
                        params.gamma1 = DynOpt.X(6);
                    end
                    
                    x_propagate = DynOpt.X;
                    for j =1:DynOpt.WindowSamples-1     %in absolute time as: k-(w-1)*Nts+1:k,
                        x_propagate =  PlantJumpMap_general_notime(x_propagate,DynOpt.model,1);
                        DynOpt.OptXstory(:,DynOpt.BackTimeIndex+j) = x_propagate;
                    end
                
                % backward optimization
                else 
                    
                    DynOpt.BackTimeIndex = k; 
                    DynOpt.temp_x0 = DynOpt.OptXstory(:,DynOpt.BackTimeIndex);
                    [NewXopt, J] = fminunc(@EvaluateCostFunctionOnWindow_general_notime_nobuf,DynOpt.temp_x0,myoptioptions);
                    
                    DynOpt.X = NewXopt;
                    DynOpt.OptXstory(:,k) = DynOpt.X;
                    
                    % params update
                    if strcmp(model,'cubli')
                        params.M = DynOpt.X(5);
                        params.If = DynOpt.X(6);
                    elseif strcmp(model,'tokamak')
                        params.gamma0 = DynOpt.X(5);
                        params.gamma1 = DynOpt.X(6);
                    end
                    
                    x_propagate = DynOpt.X;
                    for j =1:DynOpt.WindowSamples-1     %in absolute time as: k-(w-1)*Nts+1:k,
                        x_propagate =  PlantJumpMap_general_notime(x_propagate,DynOpt.model,-1);
                        DynOpt.OptXstory(:,DynOpt.BackTimeIndex-j) = x_propagate;
                    end
                end
                
                DynOpt.Jstory(k) = J;
                DynOpt.OptXstory(:,k) = DynOpt.X;
                DynOpt.OptErrorStory(:,k) = DynOpt.Xtrue - DynOpt.X;
            end
            clc;
        end
    end
    toc
    
    % storage management
    DynOpt.Y_story = DynOpt.Y_story(:,2:end);
    DynOpt.Y_full_story = DynOpt.Y_full_story(:,2:end);
    DynOpt.lambda_min = DynOpt.lambda_min(2:end);
end

