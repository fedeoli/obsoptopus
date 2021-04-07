function T_horizon = time_horizon_check(i) 

global params

if t > 1.5*T_horizon % && flag_change == 1
    T_horizon = 1.5*T;
    params.DT_burn = T_horizon/10;                                  % The impulsive DT will be applied every DT_burn. Suggestion: select odd fractions of the period
    params.N_burn_horizon = floor(T_horizon/params.DT_burn);        % Number of burns
    params.toll = [1e-3; 1e-3; 1e-3; 1e-6; 1e-6; 1e-6];
end

% Control allocation inside "params" structure
params.u = u;
params.tau = tau_PD;
params.t = time(i);

end