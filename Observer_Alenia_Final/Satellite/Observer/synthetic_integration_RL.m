%% model integration - synthetic model
function [DynOpt_out, params_out] = synthetic_integration_RL

%%%% global vars %%%%
global DynOpt params

if(DynOpt.simulationModel == 1)
    
    disp('Model simulation')
    
    % define input
    DynOpt.U = zeros(3,length(DynOpt.time));
    DynOpt.input_true = zeros(3,length(DynOpt.time));
    
    % define reference quaternion
    DynOpt.quat_ref = zeros(4,length(DynOpt.time));
    
    % free to upsate the input
    DynOpt.recollect_input = 0;
    
    % set synthetic data flag
    DynOpt.synthetic_int = 1;
    
    % store quaternions in euler angles after simulation
    DynOpt.True_quat = zeros(3,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     INTEGRATION LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init procedure
    %%% init state %%%
    if DynOpt.input_tuning == 0
        x_start = zeros(DynOpt.StateDim+length(DynOpt.param_estimate),DynOpt.Niter);
        %%% init procedure %%%%
        if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
            x_start(:,1) = [satellites_iner_ECI; DynOpt.param_estimate'];
        elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
            x_start(:,1) = [satellites_attitude; DynOpt.param_estimate'];
        elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
            x_start(:,1) = [satellites_iner_ECI; satellites_attitude; DynOpt.param_estimate'];
        end 
    else
        x_start = zeros(DynOpt.StateDim+length(DynOpt.param_estimate)+3,DynOpt.Niter);
        %%% init procedure %%%%
        if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
            x_start(:,1) = [satellites_iner_ECI; DynOpt.param_estimate';DynOpt.u_amp; DynOpt.d; DynOpt.u_freq];
        elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
            x_start(:,1) = [satellites_attitude; DynOpt.param_estimate';DynOpt.u_amp; DynOpt.d; DynOpt.u_freq];
        elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
            x_start(:,1) = [satellites_iner_ECI; satellites_attitude; DynOpt.param_estimate';DynOpt.u_amp; DynOpt.d; DynOpt.u_freq];
        end
    end
    
    params.SatellitesCoordinates = satellites_iner_ECI;
    params.SatellitesAttitude = satellites_attitude;
    
    % start integration
    tic
    disp('STARTING INTEGRATION PROCEDURE')
    for i = 2:length(time)
        
        % integration     
        [x_start(:,i), params] = model_propagate_local(i,time_step,x_start(:,i-1),params);   
        
        %%% FAULT SIMULATION %%%%
        if DynOpt.fault_sim == 1
            if i > floor(0.5*DynOpt.Niter)
                DynOpt.param_estimate = max(0.4*DynOpt.param_estimate_init,0.5*DynOpt.param_estimate);
                params_update(x_start(:,i-1));
                x_start(DynOpt.StateDim+1:end,i) = DynOpt.param_estimate;
            end
        end
        
%         DynOpt.dx_true(:,i) = AttitudeDynamics_bias_v2(x_start(:,i),params);

    end
    
    % store input
    DynOpt.input_true = DynOpt.U;
    DynOpt.desatt_true = DynOpt.desired_attitude(:,2:end);
    
    % restore initial params
    if DynOpt.fault_sim == 1
        DynOpt.param_estimate = DynOpt.param_estimate_init;
    end
    
    % state storage
    DynOpt.stateStory = zeros(DynOpt.StateDim,DynOpt.Niter);
    DynOpt.param_story = zeros(length(DynOpt.param_estimate),DynOpt.Niter);
    DynOpt.stateStory(:,1) = DynOpt.init_state;
    DynOpt.init_pos = satellites_iner_ECI;
    DynOpt.init_att = satellites_attitude;
    
    % assign storage
    DynOpt.param_story = x_start(DynOpt.StateDim+1:end,:);
    DynOpt.stateStory(:,2:end) = x_start(1:DynOpt.StateDim,2:end);
    if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
        DynOpt.position_state = DynOpt.stateStory(1:params.Nagents*6,:);
    elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
        DynOpt.attitude_state = DynOpt.stateStory(1:params.Nagents*7,:);
    elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
        DynOpt.position_state = DynOpt.stateStory(1:params.Nagents*6,:);
        DynOpt.attitude_state = DynOpt.stateStory(params.Nagents*6+1:end,:);
    end    
    
    % store quaternions in euler angles
    DynOpt.True_quat = unwrap(RotationConversion_V2_1('QtoEA321', DynOpt.attitude_state(1:4,:)')*pi/180)';
    
    %%%%%%% STORE MEASUREMENTS %%%%%
    DynOpt.eps_noise_story = DynOpt.measure_amp*randn(length(params.observed_state),DynOpt.Niter);
    
    % set synthetic data flag
    DynOpt.synthetic_int = 0;
end