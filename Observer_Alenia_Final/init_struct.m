%% simulation data init

% dependencies
addpath(genpath([pwd '/Lib']));
addpath(genpath([pwd '/Satellite']));

% clear
if exist('RL','var')
    keep RL s current_reward test_flag
    setup.load_mem = 0;
else
    clear
end
clc

% dependencies
addpath(genpath([pwd '/Lib']));
addpath(genpath([pwd '/Satellite']));

% simulation time
setup.t_start = 0;
setup.Tend = 300;
setup.Ts = 1e0;

% measurement dimension
setup.nMagneto = 2;
setup.RPYbetweenMagnetometers = 1*[0,0,90]*pi/180;
setup.dim_out = 9;

% plot and graphics
setup.plot = 0;
setup.print = 0;
setup.RL = 0;

% montecarlo or single simulation
setup.montecarlo = 0;
% set the bounds for the model generation
if setup.montecarlo == 1
    setup.lb_init = -1e-1*[1 1 1 1 1 1 1 1 1 1];
    setup.ub_init = 1e0*[1 1 1 1 1 1 1 1 1 1];
end

%%% integration %%%
setup.integration_pos = 1;
setup.integration_att = 1;

%%%%% OBSERVER %%%%%
setup.generate_plant = 1;
setup.generate_orbit = 0;
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

%%%%% OBSERVING INPUT %%%%%
setup.control = 1;
setup.u_amp = 1*pi/4;
setup.d = 0.1;
setup.T = 100;
setup.u_freq = 1/setup.T;
setup.target_attitude = 0*[1; 1; 1];
setup.lowpass_pwm = 0;

%%% FOR THE INTEGRATION TIME SEE SCNARIO_ORBITS_K.m %%%

%%%%% SAMPLING %%%%%
setup.w = 5;
setup.Nts = 3;
setup.theta = 0;
setup.beta = 0;
setup.gamma = 1;
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
setup.dJ_2 = setup.adaptive*5e-2;
setup.dJ_1 = setup.adaptive*1e-2;

%%%%% SCALE FACTOR %%%%%
setup.y_end = 3;
setup.nJ_nl = 2;
%%%% bias values %%%%%
setup.scale_factor_init = 1e0.*[1;1e-1;0*5e-2;0;5e-1].*ones(setup.y_end+setup.nJ_nl,setup.dim_out);
% memory factor
% setup.lambda = [0.8; 1; 0.7; 1; 1];
setup.lambda = 1*[1; 1; 1; 1; 1];
setup.y_weight = [1*ones(1,3) 1*ones(1,3*setup.nMagneto)];
setup.scale_factor = zeros(setup.w,setup.y_end+setup.nJ_nl,setup.dim_out);
for z = 1:setup.w
    setup.scale_factor(setup.w+1-z,:,:) = setup.scale_factor_init.*setup.lambda.^(z-1);
end
for z = 1:setup.dim_out
    setup.scale_factor(:,:,z) = setup.scale_factor(:,:,z)*setup.y_weight(z);
end
clear z

% optimisation
setup.fcon_flag = 0;
if setup.fcon_flag == 1
    setup.fmin = @fmincon;
else
    setup.fmin = @fminsearch;
end
% state bounds
setup.lb = 1*[-Inf -Inf -Inf -Inf -Inf -Inf -Inf 0];
setup.ub = 1*[Inf Inf Inf Inf Inf Inf Inf Inf];

% global search solutions
setup.globalsearch = 0;
setup.multistart = 0;
setup.nstart = 5;

%%%%% GRADIENT DESCENT %%%%%
setup.alpha_grad = 1e-11;
setup.grad_thresh = 1e-5;
setup.alpha_dyn = 1;

%%%%% INTEGRATION %%%%%
setup.forward = 1;

%%%%% BIAS %%%%%
setup.bias_dyn = 0;
setup.bias_enable = 0;
setup.bias_mag_enable = 0;
setup.optimise_params = 1;
setup.nbias = 0;
setup.nparams = 3;
setup.inertia = 1;

%%%%% NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setup.noise_enable = 1;
setup.ErrorOnEA = 1;
setup.rand_init = 0;
% bias
if setup.bias_enable 
    setup.bias = 1*ones(3,1)*5*pi/180; 
else
    setup.bias = zeros(3,1);
end
if setup.bias_mag_enable 
    setup.bias_mag = 1*abs(1e-4*randn(3*setup.nMagneto,1)+1e-3);
else
    setup.bias_mag = zeros(3*setup.nMagneto,1);
end
setup.bias_tot = setup.bias;

% measures
setup.EulerAngleNoiseOnMag = 0*setup.noise_enable*1e-2;
setup.noise_amp = 1*setup.noise_enable*[1*1e-4*ones(3,1); 1*1e-4*ones(3*setup.nMagneto,1)];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% MODEL %%%%%
setup.model = 'satellite';