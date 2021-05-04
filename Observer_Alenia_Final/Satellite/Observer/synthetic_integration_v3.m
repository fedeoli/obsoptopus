%% model integration - synthetic model
function [DynOpt,params] = synthetic_integration_v3(DynOpt,params,struct,satellites_iner_ECI,satellites_attitude)

    if(DynOpt.simulationModel == 1)

        if struct.print
            disp('Model simulation')
        end

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
            x_start(:,1) = [satellites_iner_ECI; DynOpt.param_estimate'];
        elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
            x_start(:,1) = [satellites_attitude; DynOpt.param_estimate];
        elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
            x_start(:,1) = [satellites_iner_ECI; satellites_attitude; DynOpt.param_estimate];
        end 

        params.SatellitesCoordinates = satellites_iner_ECI;
        params.SatellitesAttitude = satellites_attitude;
        params_true = params;

        %%% init magnetic field %%%
        myutc = [2019 12 15 10 20 36]; 
        LatLongAlt = eci2lla(transpose(params.SatellitesCoordinates(1:3)*1E3),myutc); 
        [mag_field_vector,~,~,~,~] = igrfmagm(max(1000,min(LatLongAlt(3),6E5)),LatLongAlt(1),LatLongAlt(2),decyear(2019,12,15),13); 
        DynOpt.mag_field_story(:,1) = mag_field_vector;

        % start integration
        if struct.print
            disp('STARTING INTEGRATION PROCEDURE')
        end
        for i = 2:length(DynOpt.time)

            % integration     
            [x_start(:,i), params_true, DynOpt] = DynOpt.model_propagate(DynOpt,i,params.time_step,x_start(:,i-1),params_true);   

        end

        % store input
        DynOpt.input_true = DynOpt.U;
        DynOpt.desatt_true = DynOpt.desired_attitude(:,2:end);

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
        DynOpt.True_quat = DynOpt.wrap(quat2eul(DynOpt.attitude_state(1:4,:)'))';

        %%%%%%% STORE MEASUREMENTS %%%%%
    %     DynOpt.eps_noise_story = DynOpt.measure_amp*randn(length(params.observed_state),DynOpt.Niter);
        DynOpt.eps_noise_story = DynOpt.measure_amp.*randn(DynOpt.dim_out,DynOpt.Niter);

        % set synthetic data flag
        DynOpt.synthetic_int = 0;
    end
end