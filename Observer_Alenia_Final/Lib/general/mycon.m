%% nonlinear constraints for fmincon
function [c, ceq] = mycon(x,DynOpt,params)
    
    % don't want nan results
    try
        n_iter = sum(DynOpt.Y_space);
        x_propagate = [params.SatellitesCoordinates; x];
        for i = 1:n_iter
            [x_propagate, params] = DynOpt.model_propagate(DynOpt,DynOpt.Niter,DynOpt.Ts,x_propagate,params);
        end
        nanFlag = isnan(x_propagate);
    catch
        nanFlag = 1;
    end
    
    
    % remember that c is forced to be <= 0
    ceq = sum(nanFlag);
    c = -1;
end