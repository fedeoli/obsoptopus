%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    ATTITUDE BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DynOpt params

if ~exist('DynOpt.ActualTimeIndex', 'var')
    iter = i;
else
    iter = DynOpt.ActualTimeIndex;
end

if params.Control
    if params.Attitude

        %%%%% Reference Attitude %%%%%
        params = ComputeReferenceAttitude_V2_4(t, params);
        DesiredAttitude_out(:,i+1) = reshape(params.DesiredAttitude, [3*(N_deputy + 1), 1]);

        %%%%% Attitude Control %%%%%
        tau_PD = AttitudeControl_V2_4(temp_att(:,end), temp_pos(:,end), t, params);
        tau_PD_out(:,i,:) = reshape(tau_PD*1e6, [3, 1, N_deputy + 1]);

        %%%%% Environmental Torques %%%%%
        [tau_GG, tau_D] = EnvironmentalTorques_V2_2(temp_att(:,end), params);
        tau_GG_out(:,i,:) = reshape(tau_GG, [3, 1, N_deputy + 1]);
        tau_D_out(:,i,:) = reshape(tau_D, [3, 1, N_deputy + 1]);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                         ACTUATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if params.Attitude
        [deputy_rel_LVLH_estimated, P] = EstimatedState_V2_1(deputy_rel_LVLH);
        [u, temp_pos(:,end), params] = SingleAxisThruster_V2_5(iter, u, temp_pos(:,end), deputy_rel_LVLH_estimated, params);
        params = CollisionAvoidanceStatusCheck_V2_1(params);
    end
end