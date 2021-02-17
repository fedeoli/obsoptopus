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
        if DynOpt.input_tuning == 0
        % inertia
%             params.sat(1).I(1,1) = state(end-2);
%             params.sat(1).I(2,2) = state(end-1);
%             params.sat(1).I(3,3) = state(end);
%             params.params_estimate = [params.sat(1).I(1,1), params.sat(1).I(2,2), params.sat(1).I(3,3)];
            
        % bias
              params.bias = state(end);
              params.params_estimate = params.bias;
        else
        % inertia
%             params.sat(1).I(1,1) = state(end-5);
%             params.sat(1).I(2,2) = state(end-4);
%             params.sat(1).I(3,3) = state(end-3);
%             params.params_estimate = [params.sat(1).I(1,1), params.sat(1).I(2,2), params.sat(1).I(3,3)];

        % bias
              params.bias = state(end-3);
              params.params_estimate = params.bias;
        end
    end
    
    params_out = params;
    
end