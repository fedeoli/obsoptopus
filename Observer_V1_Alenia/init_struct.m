%% simulation data init

% clear
if exist('i','var') && exist('n_sim','var')
    keep i n_sim
else
    clear
end
clc

% dependencies
addpath(genpath([pwd '/Lib']));
addpath(genpath([pwd '/Satellite']));

% measurement dimension
struct.dim_out = 3;

% montecarlo or single simulation
struct.montecarlo = 0;
% set the bounds for the model generation
if struct.montecarlo == 1
    struct.lb_init = -1e-1*[1 1 1 1 1 1 1 1 1 1];
    struct.ub_init = 1e0*[1 1 1 1 1 1 1 1 1 1];
end

%%%%% OBSERVER %%%%%
struct.ObserverOn = 0;
struct.OptimisationOn = 1;
struct.simulationModel = 1;
struct.identify = 1;
struct.check = 0;
struct.fault_sim = 0;
struct.always_opt = 0;
struct.flush_buffer = 1;
struct.control = 0;
struct.optimise_input = 0;
struct.optimise_params = 1;
struct.input_tuning = 0;

%%%%% OBSERVING INPUT %%%%%
struct.u_amp = 0*4e-1;
struct.d = 0.5;
struct.u_freq = 1e-6;
struct.target_attitude = 1e0*[0.5; 0.5; 0.5];

%%% FOR THE INTEGRATION TIME SEE SCNARIO_ORBITS_K.m %%%

%%%%% SAMPLING %%%%%
struct.w = 10;
struct.Nts = 3;
struct.dJcond_thresh = 5e-2;                                                  % the optimisation is run if COND > THRESH (set to 0 to set a nocare condition)
struct.theta = 0;
struct.beta = 0;
struct.gamma = 1;
% conditions to consider the optimisation result
struct.Jdot_thresh = 9e-1;
struct.blue_flag = 0;
% built in/gradient optimisation conditions
struct.J_thresh = [1e-10, 1e3];
struct.max_iter = 30;
struct.safety_density = 4;

%%%%% SCALE FACTOR %%%%%
struct.y_end = 3;
struct.nJ_nl = 1;
struct.scale_factor = 1e0.*[1;0;5e-1;2e-3].*ones(struct.y_end+struct.nJ_nl,struct.dim_out);

% optimisation
struct.fcon_flag = 0;
if struct.fcon_flag == 1
    struct.fmin = @fmincon;
else
    struct.fmin = @fminsearch;
end
% state bounds
struct.lb = 1*[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
struct.ub = Inf*[1 1 1 1 1 1 1 1 1 1];

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

%%%%% NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct.noise_enable = 0;
struct.nparams = 1;
% measures5
struct.noise_amp = 1*struct.noise_enable*1e-5;
%%%% STATE %%%%
% POS + ATT
struct.init_error_amp = 1*struct.noise_enable*[0*ones(3,1); 0*ones(3,1); 5e-2*ones(4,1); 5e-2*ones(3,1)];
% POS
% struct.init_error_amp = 0*struct.noise_enable*[1e1*ones(3,1); 5e-1*ones(3,1)];
% ATT (errors in percentage)
% struct.init_error_amp = 1*struct.noise_enable*[0; ones(3,1); ones(3,1)]*10;
%%%% PARAMS %%%%
struct.ParamsAllPos = 1;
struct.init_param_error_amp = 1*struct.noise_enable*ones(1,struct.nparams)*30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% MODEL %%%%%
struct.model = 'satellite';