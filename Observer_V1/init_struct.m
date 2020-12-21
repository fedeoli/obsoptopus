%% simulation data init

% clear
clc

% dependencies
addpath(genpath([pwd '/Lib']));

%%%%% OBSERVER %%%%%
struct.ObserverOn = 1;
struct.simulationModel = 0;
struct.identify = 0;
struct.check = 0;

% measure
struct.dim_out = 1;
struct.y_end = 1;

% optimisation
struct.opt_method = 'default';
struct.opt_algorithm = 'default';
struct.fcon_flag = 0;

%%%%% SAMPLING %%%%%
struct.w = 20;
struct.Nts = 10;

%%%%% NÌINTEGRATION %%%%%
if struct.simulationModel == 1
    struct.T0 = 0;
    struct.Ts = 5e-4;
    struct.Tend = 0.4;
else
    struct.T0 = 0.2108;
    struct.Ts = 1e-6;
    struct.Tend = 0.2110;
end

struct.int_flag = 'default';
struct.forward = 1;

%%%%% DATI %%%%%
if struct.simulationModel == 0
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

%%%%% GRADIENT DESCENT %%%%%
% struct.alpha_grad = [1e-4; 1e-5; 1e-4; 1e-4];
struct.alpha_grad = 1e-6*[1; 1; 1; 1];
struct.max_iter = 20;
%%% set 10% J buono 
struct.grad_thresh = 1e-6;

%%%%% NOISE %%%%%
struct.noise_enable = 1;
struct.noise_amp = 1*struct.noise_enable*1e-3;
struct.init_error_amp = 1*struct.noise_enable*[1e-1;5e-1];
struct.init_param_error_amp = 1*struct.noise_enable*1e-1;

%%%%% MODEL %%%%%
struct.model = 'runaway';

%%%%% SCALE FACTOR %%%%%
struct.scale_factor = [1;0].*ones(2,struct.dim_out);