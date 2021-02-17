%% model propagation
function [x_propagate, params_out] = model_propagate_local(time_instant,time_step,x_start,params)

    global DynOpt
    
    % set current_pos
    DynOpt.current_pos = time_instant;

    % integrazione estimation and real
    tspan = DynOpt.time(time_instant):time_step:DynOpt.time(time_instant) + time_step;
    
    % set initial condition if bias flag
%     if (DynOpt.bias_model == 1) && (DynOpt.synthetic_int == 0)
%        x_start(8) = x_start(end);
%        x_start(5:7) = DynOpt.Y_full_story(:,DynOpt.current_pos);
%     end

    if DynOpt.input_tuning == 0
        if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
            temp_att = params.SatellitesAttitude;
            params = DynOpt.set_input(time_instant,temp_att,x_start(1:6),params);
            temp_pos = rk4_V1_1(DynOpt.model, tspan, x_start(1:6), params);
            x_propagate = [temp_pos(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
        elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
            temp_pos = params.SatellitesCoordinates;
            params = DynOpt.set_input(time_instant,x_start(1:7),temp_pos,params);
            if DynOpt.bias_dyn == 0
                temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(1:7), params);
                x_propagate = [temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
            else
                temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(1:8), params);
                x_propagate = temp_att(:,end);
            end            
        elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
            params = DynOpt.set_input(time_instant,x_start(1:7),x_start(1:6),params);
            temp_pos = rk4_V1_1(DynOpt.model_inertial, tspan, x_start(1:6), params);
            if DynOpt.bias_dyn == 0
                temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(7:13), params);
                x_propagate = [temp_pos(:,end); temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
            else
                temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(7:14), params);
                x_propagate = [temp_pos(:,end); temp_att(:,end)];
            end  
            
        end
    else
        if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
            temp_att = params.SatellitesAttitude;
            params = DynOpt.set_input(time_instant,temp_att,x_start(1:6),params);
            temp_pos = rk4_V1_1(DynOpt.model, tspan, x_start(1:6), params);
            x_propagate = [temp_pos(:,end); x_start(end-DynOpt.nparams-2:end)];
        elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
            temp_pos = params.SatellitesCoordinates;
            params = DynOpt.set_input(time_instant,x_start(1:7),temp_pos,params);
            temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(1:7), params);
            x_propagate = [temp_att(:,end); x_start(end-DynOpt.nparams-2:end)];
        elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
            params = DynOpt.set_input(time_instant,x_start(1:7),x_start(1:6),params);
            temp_pos = rk4_V1_1(DynOpt.model(1), tspan, x_start(1:6), params);
            temp_att = rk4_V1_1(DynOpt.model(2), tspan, x_start(7:13), params);
            x_propagate = [temp_pos(:,end); temp_att(:,end); x_start(end-DynOpt.nparams-2:end)];
        end
    end
    
    params.SatellitesCoordinates = temp_pos(:,end);
    params.SatellitesAttitude = temp_att(:,end);
    
    params_out = params;
end