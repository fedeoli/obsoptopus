%% parameter update after optimisation
function [params_out,DynOpt] = params_update_local_function(state,params,DynOpt)
    
    offset = DynOpt.integration_pos*6;
    
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
        if DynOpt.bias_enable && DynOpt.optimise_params
            params.bias = state(offset+8:offset+8+DynOpt.nbias-1);
        else
            params.bias = [];
        end
        
        % inertia
        if DynOpt.inertia && DynOpt.optimise_params
  
            % 1 inertia
%             params.sat(1).I(1,1) = state(end);
%             params.inertia = params.sat(1).I(1,1);
            
            % 3 inertias
            params.sat(1).I(2,2) = state(end-1);
            params.sat(1).I(3,3) = state(end);
            params.inertia = [params.sat(1).I(1,1); params.sat(1).I(2,2); params.sat(1).I(3,3)];
        else
            params.inertia = [];
        end
            
        params.params_estimate = [params.bias; params.inertia];
    end
    
    params_out = params;
    
end