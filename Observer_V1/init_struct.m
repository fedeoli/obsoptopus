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

% montecarlo or single simulation
struct.montecarlo = 0;
if struct.montecarlo == 1
    struct.lb_init = -1e-1*[1 1 1 1];
    struct.ub_init = 1e0*[1 1 3 5];
end

%%%%% OBSERVER %%%%%
struct.ObserverOn = 0;
struct.OptimisationOn = 1;
struct.simulationModel = 1;
struct.identify = 1;
struct.check = 1;
struct.always_opt = 0;

%%%%% INTEGRATION %%%%%
if struct.simulationModel == 1
    struct.T0 = 0.2110;
    struct.Ts = 1e-6;
    struct.Tend = 0.2150;
else
    struct.T0 = 0.2110;
    struct.Ts = 1e-6;
    struct.Tend = 0.2112;
end

% measure
struct.dim_out = 1;

%%%%% SAMPLING %%%%%
struct.w = 18;
struct.Nts = 3;
struct.dJcond_thresh = 1e0;                                                % the optimisation is run if +6COND > THRESH (set to 0 to set a nocare condition)
struct.theta = 1;
% conditions to consider the optimisation result
struct.Jdot_thresh = 4e-1;
struct.blue_flag = 0;
% built in/gradient optimisation conditions
struct.J_thresh = [1e-8, 1e3];
struct.max_iter = 20;
struct.safety_density = 2;

% optimisation
struct.fcon_flag = 1;
if struct.fcon_flag == 1
    struct.fmin = @fmincon;
else
    struct.fmin = @fminunc;
end

% global search solutions
struct.globalsearch = 0;
struct.multistart = 0;
struct.nstart = 5;

%%%%% SCALE FACTOR %%%%%
struct.y_end = 3;
struct.scale_factor = 1e0.*[1;1e-11;1].*ones(struct.y_end,struct.dim_out);

%%%%% GRADIENT DESCENT %%%%%
struct.alpha_grad = 1e-11;
struct.grad_thresh = 1e-5;
struct.alpha_dyn = 1;

struct.int_flag = 'default';
struct.forward = 1;

%%%%% DATI %%%%%
if struct.check == 1
    load simulation/shots/DatiW_Shot43649.mat
    %%% filtro %%%
    struct.filter = 0;
    if struct.filter == 1
        temp_dati = moving_average_v2(dati(:,2),15,'center');
        struct.dati = [dati(:,1),temp_dati];
        struct.dati_raw = dati;
    else
        struct.dati = dati;
        struct.dati_raw = dati;
    end
else
    struct.filter = 0;
end
%%%%%%%%%%%%%%

% sample time
struct.sample_time = 1;

%%%%% NOISE %%%%%
struct.noise_enable = 0;
struct.ParamsAllPos = 1;
struct.noise_amp = 1*struct.noise_enable*5e-4;
struct.init_error_amp = 1*struct.noise_enable*[20;20;0;0];
struct.init_param_error_amp = 1*struct.noise_enable*[80;80];

%%%%% MODEL %%%%%
struct.model = 'runaway';

% state bounds
struct.lb = [-eps -eps -Inf -Inf 0 0 0 0];
struct.ub = Inf*[1 1 1 1 1 1 1 1];