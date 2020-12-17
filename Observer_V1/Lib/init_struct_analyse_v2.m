%% simulation data init
function struct = init_struct_analyse_v2(w,Nts)
    %% simulation data init

    %%%%% OBSERVER %%%%%
    struct.ObserverOn = 1;
    struct.opt_method ='default';
    struct.simulationModel = 0;
    struct.y_end = 2;

    %%%%% SAMPLING %%%%%
    struct.w = w;
    struct.Nts = Nts;

    %%%%% NÌINTEGRATION %%%%%
    if struct.simulationModel == 1
        struct.T0 = 0;
        struct.Ts = 5e-4;
        struct.Tend = 1.5;
    else
        struct.T0 = 0.2108;
        struct.Ts = 5e-4;
        struct.Tend = 0.2130;
    end

    struct.forward = 1;

    %%%%% DATI %%%%%
    load simulation/shots/DatiW_Shot43649.mat
    struct.dati = dati;
    struct.sample_time = 1;

    %%%%% GRADIENT DESCENT %%%%%
    % struct.alpha_grad = [1e-4; 1e-5; 1e-4; 1e-4];
    struct.alpha_grad = 1e-6*[1; 1; 1; 1];
    struct.max_iter = 100;
    %%% set 10% J buono 
    struct.grad_thresh = 1e-6;

    %%%%% NOISE %%%%%
    struct.noise_enable = 1;
    struct.noise_amp = struct.noise_enable*1e-4;
    struct.init_error_amp = struct.noise_enable*[5e-2; 5e-2];
    struct.init_param_error_amp = struct.noise_enable*5e-1;

    %%%%% MODEL %%%%%
    struct.model = 'runaway';

    %%%%% PESI Y %%%%%
    % struct.weight = 1*ones(1,struct.w);
    struct.weight = 1*linspace(0.1,1,struct.w);
    % struct.weight = 1*linspace(1,0.11,struct.w);

    %%%%% PESI dY %%%%%
    struct.dweight = 0*ones(1,struct.w);
    % struct.dweight = 1e-1*linspace(1,0.1,struct.w);
    % struct.dweight = 1e-1*linspace(0.1,1,struct.w);

    %%%%% PESI ddY %%%%%
    struct.ddweight = 0*ones(1,struct.w);

    %%%%% PESI intY %%%%%
    struct.intweight = 0*ones(1,struct.w);

    %%%%% SCALE FACTOR %%%%%
    struct.scale_factor = [1 1e-2 1 1];
end