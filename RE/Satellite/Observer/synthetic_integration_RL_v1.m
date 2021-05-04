%% model integration - synthetic model
function [DynOpt_out,params_out] = synthetic_integration_RL_v1

    %%% global vars %%%
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
        x_start = zeros(DynOpt.StateDim+length(DynOpt.param_estimate),DynOpt.Niter);
        %%% init procedure %%%%
        if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
            x_start(:,1) = [params.SatellitesCoordinates; DynOpt.param_estimate'];
        elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
            x_start(:,1) = [params.SatellitesAttitude; DynOpt.param_estimate'];
        elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
            x_start(:,1) = [params.SatellitesCoordinates; params.SatellitesAttitude; DynOpt.param_estimate'];
        end
        
        % init output
        Y_true = DynOpt.get_measure(x_start(:,1),params);

        % start integration
        tic
        disp('STARTING INTEGRATION PROCEDURE')
        for i = 2:DynOpt.Niter    
            % Display iteration step
            if mod(i,10) == -1
                clc
                disp(['Iteration Number: ', num2str(i),'/',num2str(DynOpt.Niter-1)])
            end
            
            % integration     
            [x_start(:,i), params] = model_propagate_local_RL_v1(i,DynOpt.tspan,x_start(:,i-1),params);   
            Y_true(:,i) = DynOpt.get_measure(x_start(:,i),params);
        end

        % store input
        DynOpt.input_true = DynOpt.U;
        DynOpt.desatt_true = DynOpt.desired_attitude(:,2:end);

        % state storage
        DynOpt.stateStory = zeros(DynOpt.StateDim,DynOpt.Niter);
        DynOpt.param_story = zeros(length(DynOpt.param_estimate),DynOpt.Niter);
        DynOpt.stateStory(:,1) = DynOpt.init_state;
        DynOpt.init_pos = params.SatellitesCoordinates;
        DynOpt.init_att = params.SatellitesAttitude;

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
        
        % store output
        DynOpt.Ytrue = Y_true;

        % store quaternions in euler angles
        DynOpt.True_quat = DynOpt.wrap(quat2eul(DynOpt.attitude_state(1:4,:)'))';

        %%%%%%% STORE MEASUREMENTS %%%%%
        DynOpt.eps_noise_story = DynOpt.measure_amp*randn(DynOpt.dim_out,DynOpt.Niter);

        % set synthetic data flag
        DynOpt.synthetic_int = 0;
    end
    
    % output
    DynOpt_out = DynOpt;
    params_out = params;
    
end