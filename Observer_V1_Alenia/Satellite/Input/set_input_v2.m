%% Set input function
function params_out = set_input_v2(pos,x_att,x_pos,params)

    global DynOpt
    
    % system control input
    if DynOpt.params.input_flag == 1
        %%% position input
        params.u = 0;
        
        %%% attitude input
        n_act = length(DynOpt.target_attitude);
        t = DynOpt.time(pos);
        real_module = floor(t./(2*pi*DynOpt.u_freq));
        switch_pwm = 1-2*mod(real_module,2);
        temp_att = DynOpt.target_attitude + DynOpt.u_amp.*switch_pwm.*ones(n_act,1);
        params.DesiredAttitude =  temp_att;
    
        % compute input
        tau_PD = AttitudeControl_V2_5(x_att, x_pos, DynOpt.time(pos), params);
        params.tau = tau_PD;
    else
        % position input
        params.u = 0;
       
        %%% attitude input
        params.tau = zeros(3,1);
    end
    
    params_out = params;
end