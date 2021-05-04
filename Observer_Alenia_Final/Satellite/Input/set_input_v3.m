%% Set input function
function params_out = set_input_v3(pos,x_att,x_pos,params)

    global DynOpt
    
    % system control input
    if DynOpt.params.input_flag == 1
        %%% position input
        params.u = 0;
        
        %%% attitude input
        n_act = length(DynOpt.target_attitude);
        t = DynOpt.time(pos);
        
        % set: T
        T = floor(DynOpt.Ts./(DynOpt.d)./DynOpt.u_freq);
        real_module = floor(t./T);
        % duty cycle (ceil to avoid k=1???)
        k = floor((DynOpt.d).*T);
        % module
        mod_u = mod(real_module,T);
        
        if mod_u < k
            switch_pwm = 1;
        else
            switch_pwm = -1;
        end
        
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