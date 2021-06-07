%% function for RL (sat)
function [setup,RL] = save_past_sat(setup,DynOpt,DynOpt_save,RL)
    %%% set correct initial condition for next step %%%
    offset = setup.integration_pos*6;
    RL.S.S0 = DynOpt_save.OptXstory_runtime(offset+1:end,end);
    RL.S.satellites_iner_ECI_true = DynOpt_save.position_state(:,end);
    RL.S.satellites_attitude_true = DynOpt_save.attitude_state(:,end);

    % disable initial error on the estmation, we get the last opt state
%     setup.init_error_amp = 0*setup.init_error_amp;
%     setup.init_param_error_amp = 0*setup.init_param_error_amp;

    setup.t_start = setup.Tend;
    setup.Tend = setup.t_start + setup.T_duration;

    % save buffer mems
    setup.load_mem = 1;
    setup.Y_last = DynOpt.Y_last;
    setup.dY_last = DynOpt.dY_last;
    setup.intY_last = DynOpt.intY_last;
    setup.buf_dY_last = DynOpt.buf_dY_last;
    setup.Y_full_story_last = DynOpt.Y_full_story_last;
    setup.dY_full_story_last = DynOpt.dY_full_story_last;
    setup.DynOpt.intY_full_story_last = DynOpt.intY_full_story_last;

    setup.buf_dYhat_last = DynOpt.buf_dYhat_last;
    setup.Yhat_full_story_last = DynOpt.Yhat_full_story_last;
    setup.dYhat_full_story_last = DynOpt.dYhat_full_story_last;
    setup.DynOpt.intYhat_full_story_last = DynOpt.intYhat_full_story_last;

    setup.Y_space_last = DynOpt.Y_space_last;
    setup.Y_space_full_story_last = DynOpt.Y_space_full_story_last;

    setup.OptXstoryTRUE_last = DynOpt.OptXstoryTRUE_last;
    setup.OptXstory_last = DynOpt.OptXstory_last;
    setup.OptXstory_runtime_last = DynOpt.OptXstory_runtime_last;
    setup.Xstory_last = DynOpt.Xstory_last;

    setup.time_last = DynOpt.time_last;

    setup.input_true_last = DynOpt.input_true_last;

    setup.mag_field_story_last = DynOpt.mag_field_story_last;

    setup.eps_noise_story_last = DynOpt.eps_noise_story_last;
    
    % track error
    setup.buf_trackerr_last = DynOpt.buf_trackerr;
end