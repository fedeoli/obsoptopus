%% parameter update after optimisation
function params_out = params_update_local(state,params)

    global DynOpt
    
    if strcmp(DynOpt.modelname,'cubli')
        params.M = state(5);
        params.If = state(6);
    elseif strcmp(DynOpt.modelname,'tokamak')
        params.gamma0 = state(5);
        params.gamma1 = state(6);
    elseif strcmp(DynOpt.modelname,'runaway')
        params.gamma = state(3);
        params.gamma1 = state(4);
    elseif strcmp(DynOpt.modelname,'satellite')      
        % bias
              params.bias = state(end);
              params.params_estimate = params.bias;
        
        % inertia
%             params.sat(1).I(1,1) = state(end-2);
%             params.sat(1).I(2,2) = state(end-1);
%             params.sat(1).I(3,3) = state(end);
%             params.params_estimate = [params.sat(1).I(1,1), params.sat(1).I(2,2), params.sat(1).I(3,3)];
    end
    
    params_out = params;
    
end