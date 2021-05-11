%% model propagation
function [x_propagate, params_out, DynOpt] = model_propagate_local_v3(DynOpt,time_instant,time_step,x_start,params)
    
    % set current_pos
    DynOpt.current_pos = time_instant;

    % integrazione estimation and real
    tspan = DynOpt.time(1):time_step:DynOpt.time(2) + time_step;

    if DynOpt.integration_pos == 1 && DynOpt.integration_att == 0
        % current attitude 
        temp_att = params.SatellitesAttitude;
        % set input
        [params,DynOpt] = DynOpt.set_input(DynOpt,time_instant,temp_att,x_start(1:6),params);
        % propagate
        [temp_pos,DynOpt] = rk4_V1_1_function(DynOpt.model, tspan, x_start(1:6), params,DynOpt);
        % wrap and update
        x_propagate = [temp_pos(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
        params.SatellitesCoordinates = temp_pos(:,end);
    elseif DynOpt.integration_pos == 0 && DynOpt.integration_att == 1
        % current orbit %
        temp_pos = params.SatellitesCoordinates;
        % set input
        [params,DynOpt] = DynOpt.set_input(DynOpt,time_instant,x_start(1:7),temp_pos,params);
        % propagate attitude
        if DynOpt.bias_dyn == 0
            [temp_att,DynOpt] = rk4_V1_1_function(DynOpt.model, tspan, x_start(1:7), params, DynOpt);
            x_propagate = [temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
        else
            [temp_att,DynOpt] = rk4_V1_1(DynOpt.model, tspan, x_start(1:end), params, DynOpt);
            x_propagate = [temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-2):end)];
        end  
        % wrap and update
        params.SatellitesAttitude = temp_att(:,end);
    elseif DynOpt.integration_pos == 1 && DynOpt.integration_att == 1
        % set input
        [params,DynOpt] = DynOpt.set_input(DynOpt,time_instant,x_start(7:13),x_start(1:6),params);
        % propagate orbit
        [temp_pos,DynOpt] = rk4_V1_1_function(DynOpt.model_inertial, tspan, x_start(1:6), params,DynOpt);
        params.SatellitesCoordinates = temp_pos(:,end);
        % propagate attitude
        if DynOpt.bias_dyn == 0
            [temp_att,DynOpt] = rk4_V1_1_function(DynOpt.model, tspan, x_start(7:13), params,DynOpt);
            x_propagate = [temp_pos(:,end); temp_att(:,end); x_start(end-(length(DynOpt.param_estimate)-1):end)];
        else
            [temp_att,DynOpt] = rk4_V1_1_function(DynOpt.model, tspan, x_start(7:end), params,DynOpt);
            x_propagate = [temp_pos(:,end); temp_att(:,end)];
        end  
        % wrap and update
        params.SatellitesAttitude = temp_att(:,end);
    end
    
    
    %%%%%%%% store magnetic field %%%%%%%
    if DynOpt.synthetic_int
        myutc = [2019 12 15 10 20 36]; %CHANGE THIS...??
        LatLongAlt = eci2lla(transpose(params.SatellitesCoordinates(1:3)*1E3),myutc); %converto from ECI to latitude, longitude,  altitude
        [mag_field_vector,~,~,~,~] = igrfmagm(max(1000,min(LatLongAlt(3),6E5)),LatLongAlt(1),LatLongAlt(2),decyear(2019,12,15),13); %mag_field_vector is in nanotesla, by IGRF11-12
        DynOpt.mag_field_story(:,time_instant) = mag_field_vector;
    end
    
    %%%% output params %%%%
    params_out = params;
end