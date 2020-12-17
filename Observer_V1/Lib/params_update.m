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
        params.gamma = state(3);
        params.gamma1 = state(4);
%         params.ni = state(5);
%         params.Wt = state(6);
%         params.Q = state(7);
%         params.S = state(8);
    elseif strcmp(DynOpt.modelname,'satellite')
        params.sat(1).M = state(end);
        params.params_estimate = params.sat(1).M;
    end
    
end