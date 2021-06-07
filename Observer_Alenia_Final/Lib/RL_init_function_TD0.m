%% update DynOpt with RL data 
function [DynOpt,satellites_iner_ECI,satellites_attitude] = RL_init_function_TD0(DynOpt,params,satellites_iner_ECI,struct)

    %%%%%%%%%%%%%%%%%%% SET INITIAL CONDITIONS FOR GREEDY %%%%%%%%%%%%%%%%%
    RL = DynOpt.RL;

    %%%%%%%%%%%%% from the actions %%%%%%%%%%%%%%%%%%%    
    %%% Magnetometers %%%
    DynOpt.nMagneto = struct.nMagneto;%*RL.S.A(1,RL.S.i); 
    DynOpt.dim_out = 3+3*DynOpt.nMagneto;
    DynOpt.measure_amp = DynOpt.measure_amp(1:DynOpt.dim_out);

    %%% Aw %%%
    temp_Aw = zeros(3,1);
    for i=1:3
        if RL.S.A(RL.S.i) == 0
            temp_Aw(i) = 0;
        else
            if RL.A.Aw_sign(i) ~= 0
                temp_Aw(i) = -RL.S.A(RL.S.i)*RL.A.Aw_sign(i);
            else
                temp_Aw(i) = RL.S.A(RL.S.i);
            end
        end
    end
    DynOpt.Aw = 0*RL.A.step*temp_Aw;

    %%%%%%%%%%%%%%% from the state %%%%%%%%%%%%%%%%%%%%%%%%
    % generate state init
    if RL.S.i > 1
        % true iner_ECI
        satellites_iner_ECI = RL.S.satellites_iner_ECI_true; 
        % true attitude
        satellites_attitude = RL.S.satellites_attitude_true;
    else
        % true attitude
        RL.S.satellites_attitude_true = RL.S.satellites_attitude_singleopt;
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

    DynOpt = scale_factor(DynOpt);
end