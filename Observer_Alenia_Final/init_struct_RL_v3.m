%% simulation data init

% clear
if exist('RL','var')
    keep RL
else
    clear
end
clc

% dependencies
addpath(genpath([pwd '/Lib']));
addpath(genpath([pwd '/Satellite']));

% simulation time
struct.Tend = 250;
struct.Ts = 1e0;

% measurement dimension
struct.nMagneto = 2;
struct.RPYbetweenMagnetometers = 1*[0,0,90]*pi/180;
struct.dim_out = 3+3*struct.nMagneto;

% plot and graphics
struct.plot = 1;
struct.print = 1;

% montecarlo or single simulation
struct.montecarlo = 0;
% set the bounds for the model generation
if struct.montecarlo == 1
    struct.lb_init = -1e-1*[1 1 1 1 1 1 1 1 1 1];
    struct.ub_init = 1e0*[1 1 1 1 1 1 1 1 1 1];
end

%%% integration %%%
struct.integration_pos = 1;
struct.integration_att = 1;

%%%%% OBSERVER %%%%%
struct.generate_plant = 1;
struct.ObserverOn = 1;
struct.Observer = 'OPT';
struct.OptimisationOn = 1;
struct.simulationModel = 1;
struct.identify = 1;
struct.check = 0;
struct.fault_sim = 0;
struct.flush_buffer = 0;
struct.always_opt = 0;
struct.control = 1;
struct.optimise_input = 0;
struct.input_tuning = 0;

%%%%% OBSERVING INPUT %%%%%
struct.u_amp = 1*pi/6;
struct.d = 0.1;
struct.T = 100;
struct.u_freq = 1/struct.T;
struct.target_attitude = pi/4*[1; 1; 1];
struct.lowpass_pwm = 0;

%%% FOR THE INTEGRATION TIME SEE SCNARIO_ORBITS_K.m %%%

%%%%% SAMPLING %%%%%
struct.w = 5;
struct.Nts = 3;
struct.theta = 0;
struct.beta = 0;
struct.gamma = 1;
% conditions to consider the optimisation result
struct.Jdot_thresh = 9e-1;
struct.blue_flag = 0;
% built in/gradient optimisation conditions
struct.J_thresh = [1e-10, 1e3];
struct.max_iter = 40;
struct.maxFcount = Inf;
struct.safety_density = 10;

%%%% HYSTERESIS %%%%
% the optimisation is run if COND > THRESH (set to 0 for a nocare condition)
struct.dJ_1 = 2e-2;
struct.dJ_2 = 5e-2;

%%%%% SCALE FACTOR %%%%%
struct.y_end = 3;
struct.nJ_nl = 2;
struct.scale_factor_init = 1e0.*[1;1e-1;1e-2;0;5e-1].*ones(struct.y_end+struct.nJ_nl,struct.dim_out);
% memory factor
% struct.lambda = [0.8; 1; 0.7; 1; 1];
struct.lambda = 1*[1; 1; 1; 1; 1];
struct.y_weight = [1*ones(1,3) 1*ones(1,3*struct.nMagneto)];
struct.scale_factor = zeros(struct.w,struct.y_end+struct.nJ_nl,struct.dim_out);
for z = 1:struct.w
    struct.scale_factor(struct.w+1-z,:,:) = struct.scale_factor_init.*struct.lambda.^(z-1);
end
for z = 1:struct.dim_out
    struct.scale_factor(:,:,z) = struct.scale_factor(:,:,z)*struct.y_weight(z);
end
clear z

% optimisation
struct.fcon_flag = 0;
if struct.fcon_flag == 1
    struct.fmin = @fmincon;
else
    struct.fmin = @fminsearch;
end
% state bounds
struct.lb = 1*[-Inf -Inf -Inf -Inf -Inf -Inf -Inf 0];
struct.ub = 1*[Inf Inf Inf Inf Inf Inf Inf Inf];

% global search solutions
struct.globalsearch = 0;
struct.multistart = 0;
struct.nstart = 5;

%%%%% GRADIENT DESCENT %%%%%
struct.alpha_grad = 1e-11;
struct.grad_thresh = 1e-5;
struct.alpha_dyn = 1;

%%%%% INTEGRATION %%%%%
struct.forward = 1;

%%%%% BIAS %%%%%
struct.bias_dyn = 1;
struct.bias_enable = 1;
struct.bias_mag_enable = 0;
struct.optimise_params = 1;
struct.nparams = 3;

%%%%% NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct.noise_enable = 1;
struct.ErrorOnEA = 1;
struct.rand_init = 0;
% bias
struct.bias = 1*struct.noise_enable*struct.bias_enable*ones(3,1)*5*pi/180; 
struct.bias_mag = 1*struct.noise_enable*struct.bias_mag_enable*abs(1e-4*randn(3*struct.nMagneto,1)+1e-3);
struct.bias_tot = struct.bias; %struct.bias_mag];
% measures
struct.EulerAngleNoiseOnMag = 0*struct.noise_enable*1e-2;
struct.noise_amp = 1*struct.noise_enable*[1*1e-4*ones(3,1); 1*1e-4*ones(3*struct.nMagneto,1)];
%%%% STATE %%%%
if struct.integration_pos == 1 && struct.integration_att == 1
    % POS + ATT
    struct.init_error_amp = 1*struct.noise_enable*[0*ones(3,1); 0*ones(3,1); 40*ones(4,1); 10*ones(3,1)];
elseif struct.integration_pos == 1 && struct.integration_att == 0
    % POS
    struct.init_error_amp = 0*struct.noise_enable*[1e1*ones(3,1); 5e-1*ones(3,1)];
else
    % ATT (errors in percentage)
    struct.init_error_amp = 1*struct.noise_enable*[20*ones(4,1); 10*ones(3,1)];
end
%%%% PARAMS %%%%
struct.ParamsAllPos = 1;
struct.init_param_error_amp = 1*struct.noise_enable*ones(1,struct.nparams)*50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% MODEL %%%%%
struct.model = 'satellite';