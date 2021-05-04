%% model propagation
function x_propagate = model_propagate(time_instant,time_step,x_start)

    global params DynOpt
    
    % set current_pos
    DynOpt.current_pos = time_instant;

    % integrazione estimation and real
    tspan = DynOpt.time(time_instant):time_step:DynOpt.time(time_instant) + time_step;

    if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
        temp_att = params.SatellitesAttitude;
        DynOpt.set_input(time_instant,temp_att,x_start(1:6));
        temp_pos = rk4_V1_1(DynOpt.model, tspan, x_start(1:6), params);
        x_propagate = [temp_pos(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
    elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
        temp_pos = params.SatellitesCoordinates;
        DynOpt.set_input(time_instant,x_start(1:7),temp_pos);
        if DynOpt.bias_dyn == 0
            temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(1:7), params);
            x_propagate = [temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
        else
            temp_att = rk4_V1_1(DynOpt.model, tspan, x_start(1:8), params);
            x_propagate = temp_att(:,end);
        end
    elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
        set_input(time_instant,x_start(1:7),x_start(1:6));
        temp_pos = rk4_V1_1(DynOpt.model(1), tspan, x_start(1:6), params);
        temp_att = rk4_V1_1(DynOpt.model(2), tspan, x_start(7:13), params);
        x_propagate = [temp_pos(:,end); temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
    end
    
    if DynOpt.input_tuning == 1
       x_propagate = [x_propagate; x_start(end-2:end)]; 
    end
    
    params.SatellitesCoordinates = temp_pos(:,end);
    params.SatellitesAttitude = temp_att(:,end);
end