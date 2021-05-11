%% update DynOpt with RL data 
function [DynOpt,satellites_iner_ECI,satellites_attitude] = RL_init_function_TD0(DynOpt,params,satellites_iner_ECI)

    %%%%%%%%%%%%%%%%%%% SET INITIAL CONDITIONS FOR GREEDY %%%%%%%%%%%%%%%%%
    RL = DynOpt.RL;

    %%%%%%%%%%%%% from the actions %%%%%%%%%%%%%%%%%%%
    %%% Magnetometers %%%
    DynOpt.nMagneto = RL.S.A(1,RL.S.i);
    DynOpt.dim_out = 3+3*DynOpt.nMagneto;
    DynOpt.measure_amp = DynOpt.measure_amp(1:DynOpt.dim_out);

    %%%%%%%%%%%%%%% from the state %%%%%%%%%%%%%%%%%%%%%%%%
    % generate state init
    DynOpt.target_attitude = RL.S.T0;
    if RL.S.i > 1
        % true iner_ECI
        satellites_iner_ECI = RL.S.satellites_iner_ECI_true; 
        % true attitude
        satellites_attitude = RL.S.satellites_attitude_true;
    else
        % true attitude
        attitude_quat = eul2quat(transpose(RL.S.satellites_attitude_true(1:3)));
        satellites_attitude = [transpose(attitude_quat); RL.S.satellites_attitude_true(4:6)];
    end
    

    %%% state data %%%%
    if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
        DynOpt.StateDim = 6*params.Nagents;
        DynOpt.StateDim_single_agent = 6;
        DynOpt.init_state = satellites_iner_ECI;
    elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
        DynOpt.StateDim = 7*params.Nagents;
        DynOpt.StateDim_single_agent = 7;
        DynOpt.init_state = satellites_attitude;
    elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
        DynOpt.StateDim = 13*params.Nagents;
        DynOpt.StateDim_single_agent = 13;
        DynOpt.init_state = [satellites_iner_ECI; satellites_attitude];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% reset scale_factor depending on the new buffers %%%%%%%%%%%%
    %%%%% SCALE FACTOR %%%%%
    DynOpt.y_end = 3;
    DynOpt.nJ_nl = 2;
    %%%% bias values %%%%%
    DynOpt.scale_factor_init = 1e0.*[1;1e-1;1*5e-2;0;5e-1].*ones(DynOpt.y_end+DynOpt.nJ_nl,DynOpt.dim_out);
    % memory factor
    % struct.lambda = [0.8; 1; 0.7; 1; 1];
    DynOpt.lambda = 1*[1; 1; 1; 1; 1];
    DynOpt.y_weight = [1*ones(1,3) 1*ones(1,3*DynOpt.nMagneto)];
    DynOpt.scale_factor = zeros(DynOpt.w,DynOpt.y_end+DynOpt.nJ_nl,DynOpt.dim_out);
    for z = 1:DynOpt.w
        DynOpt.scale_factor(DynOpt.w+1-z,:,:) = DynOpt.scale_factor_init.*DynOpt.lambda.^(z-1);
    end
    for z = 1:DynOpt.dim_out
        DynOpt.scale_factor(:,:,z) = DynOpt.scale_factor(:,:,z)*DynOpt.y_weight(z);
    end
    clear z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% set orbit %%%%%
    if DynOpt.generate_orbit
        DynOpt.orbit = RL.S.orbit(:,RL.S.i);
    end
end