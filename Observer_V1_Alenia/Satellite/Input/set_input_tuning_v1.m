%% Set input function
function params_out = set_input_tuning_v1(pos,x_att,x_pos,params,theta_u)

    global DynOpt
    
    u_amp = theta_u(1);
    u_d = theta_u(2);
    u_freq = theta_u(3);
    
    % system control input
    if DynOpt.params.input_flag == 1
        
        %%% position input
        params.u = zeros(3,1);
        
        % Check the recollect_flag to see if you can store the input
        if DynOpt.recollect_input == 0 || DynOpt.optimise_input == 1
            
            %%% attitude input
            n_act = length(DynOpt.target_attitude);
            t = DynOpt.time(pos);

            % set: T
            T = floor(DynOpt.Ts./(u_d)./u_freq);
            real_module = floor(t./T);
            % duty cycle (ceil to avoid k=1???)
            k = floor((u_d).*T);
            % module
            mod_u = mod(real_module,T);

            % store past
            switch_pwm_past = DynOpt.switch_pwm;

            % assign new one
            if mod_u < k
                DynOpt.switch_pwm = 1;
            else
                DynOpt.switch_pwm = -1;
            end

            % desired attitude
            temp_att = DynOpt.target_attitude + u_amp.*DynOpt.switch_pwm.*ones(n_act,1);
            temp_att = min(pi/2, max(-pi/2, temp_att));

            % low pass filter
            if (DynOpt.switch_pwm*switch_pwm_past < 0)
                DynOpt.t_lowpass = 0;
            end
            DynOpt.t_lowpass = min(DynOpt.Niter,DynOpt.t_lowpass+1);
            time = DynOpt.time(DynOpt.t_lowpass);
            tau = 10;
            K = 0.8*tau;
            temp_att = temp_att*(1-DynOpt.switch_pwm*(K/tau)*exp(-time/tau));
            params.DesiredAttitude =  temp_att;

            % compute input
            tau_PD = AttitudeControl_V2_5(x_att, x_pos, DynOpt.time(pos), params);
            params.tau = tau_PD;
            
            % store the input and the desired attitude
            DynOpt.desired_attitude(:,pos) = params.DesiredAttitude;
            DynOpt.U(:,pos) = tau_PD;
        else           
            %%% attitude input
            params.tau = DynOpt.input_true(:,pos);
            
            % store the input
            DynOpt.U(:,pos) = params.tau;
        end
    else
        %%% position input
        params.u = zeros(3,1);
       
        %%% attitude input
        params.tau = zeros(3,1);
        
        params.DesiredAttitude = zeros(3,1);
        DynOpt.desired_attitude(:,pos) = params.DesiredAttitude;
    end
    
    params_out = params;
end