
%% parameter update after optimisation
function params_update(state)

    global DynOpt params
    
    if strcmp(DynOpt.modelname,'cubli')
        params.M = state(5);
        params.If = state(6);
    elseif strcmp(DynOpt.modelname,'tokamak')
        params.gamma0 = state(5);
        params.gamma1 = state(6);
    elseif strcmp(DynOpt.modelname,'runaway')
        % gamma and gamma1
%         params.gamma = state(end-5);
%         params.gamma1 = state(end-4);       
%         % Q and eps_coef
        params.Q = state(end-3);
        params.eps_coef = state(end-2);
        % ringing
        params.wq = state(end-1);
        params.chi = state(end);
    elseif strcmp(DynOpt.modelname,'satellite')
        params.sat(1).M = state(end);
        params.params_estimate = params.sat(1).M;
    end
    
end