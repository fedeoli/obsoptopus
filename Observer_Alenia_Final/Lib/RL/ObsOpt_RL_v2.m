function [DynOpt, params] = ObsOpt_RL_v2(struct,RL)

    %% Init Section
    close all
%     clc

    %%%%%% global vars %%%%%%
    global DynOpt params

    % dependencies
    addpath(genpath([pwd '/Lib']));
    addpath(genpath([pwd '/Satellite']));

    % new random seed
    % RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

    %%%% INITIAL SETUP %%%%
    DynOpt.generate_plant = 1;
    DynOpt.bias_dyn = struct.bias_dyn;
    DynOpt.bias_enable = struct.bias_enable;
    DynOpt.modelname = struct.model;
    DynOpt.ObserverOn = 0;
    DynOpt.OptimisationOn = 0;
    DynOpt.simulationModel = struct.simulationModel;
    DynOpt.nMagneto = struct.nMagneto;
    DynOpt.noise_enable = struct.noise_enable;
    DynOpt.measure_amp = struct.noise_amp;

    %%% ENABLE INPUT TUNING MODE %%%
    DynOpt.dJcond_thresh = 0;

    %%% INPUT SETUP%%%
    DynOpt.set_input = @set_input_v5;
    DynOpt.switch_pwm = 1;
    DynOpt.t_lowpass = 1;
    DynOpt.lowpass_pwm = struct.lowpass_pwm;
    
    %%% GENERAL FUNCTIONS %%%%
    DynOpt.get_measure = @get_measure_RL_v1;
    DynOpt.wrap = @wrapToPi;
    
    %%%%%%%%%%%%%%%%%%% SET INITIAL CONDITIONS FOR GREEDY %%%%%%%%%%%%%%%%%
    DynOpt.u_freq = RL.S.F0(RL.S.i);
    DynOpt.target_attitude = RL.S.T0(:,RL.S.i);
    DynOpt.u_amp = RL.S.A(1,RL.S.i);
    DynOpt.d = RL.S.A(2,RL.S.i);

    %%%% satellite init %%%%
    satellite_init;

    %%% OBSERVABILITY SYMBOLIC ANALYSIS %%
    DynOpt.ObsTol = 5e-2;
    if RL.S.i == 1
        [DynOpt,params] = SymAnalysis_RL_v1;
        [DynOpt.theta,DynOpt.dtheta,DynOpt.dtheta_num] = ObsAnalysis_RL_v2(DynOpt,2,0,0,1);
    else
        DynOpt.dtheta = RL.S.DynOpt(1).dtheta;
    end

    %% PLANT model and synthetic data
    %%% DEFINE MODEL NAME %%%
    if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
        DynOpt.model = @InertialDynamicsIntegrator_V2_2;
    elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
        DynOpt.model = @AttitudeDynamics_bias_v2;
    elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
        DynOpt.model = @AttitudeDynamics_bias_v2;
        DynOpt.model_inertial = @InertialDynamicsIntegrator_V2_2;
    end
    
    %%% initial condition %%%
    params.SatellitesAttitude = [eul2quat(RL.S.S0(1:3,RL.S.i)')';RL.S.S0(4:6,RL.S.i)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLANT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DynOpt.generate_plant == 1
        tic
        [DynOpt,params] = synthetic_integration_RL_v1;
        toc
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% OBSERVABILITY NUMERICAL DATA %%%%%%%%%%%%%%%%%%
    for i=1:DynOpt.Niter
       state = DynOpt.attitude_state(:,i);
       output = DynOpt.Ytrue(:,i);
       [~,~,DynOpt.dtheta_num] = ObsAnalysis_RL_v2(DynOpt,3,state,output,0);
%        DynOpt.dtheta_num_story(i,:,:) = DynOpt.dtheta_num;
       DynOpt.dtheta_rank_story(i) = rank(DynOpt.dtheta_num*DynOpt.dtheta_num',DynOpt.ObsTol);
       DynOpt.dtheta_eig_story_min(i) = min(eig(DynOpt.dtheta_num*DynOpt.dtheta_num'));
       DynOpt.dtheta_eig_story_max(i) = max(eig(DynOpt.dtheta_num*DynOpt.dtheta_num'));
    end
end