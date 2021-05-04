%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = get_measure_v5ECI(x,j,forward,buf_dy,buf_intY,params,current_pos)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt

x_propagate = x;

offset = DynOpt.integration_pos*6;

% check if measure or propagate
if j == 0
    measure_flag = 1;
end



%%%% FIRST OBSERVATION
if measure_flag == 1
    %%%% FIRST OBSERVATION
    x_read = zeros(DynOpt.dim_out,1);
    for i=1:length(params.observed_state)
        x_read(i) = x_propagate(params.observed_state(i));
    end 
    
    % get magnetometers
    if DynOpt.nMagneto > 0
        [MagTemp_body1, MagTemp_body2,Mag_ECI] = Magnetometers_GetMeasures_v2(DynOpt,params.SatellitesCoordinates(1:3)',x_propagate(offset+1:offset+4)',DynOpt.RPYbetweenMagSensors');
    end
    
    if DynOpt.nMagneto == 1
        Mag = Mag_ECI;
    elseif DynOpt.nMagneto == 2
        Mag = [Mag_ECI; Mag_ECI];
    else
        Mag = [];
    end

    % output computation
    x_read(4:end) = Mag;
    y_read = x_read;
    
    % derivatives
    dy = zeros(DynOpt.dim_out,1);
    for k=1:DynOpt.dim_out
        [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,y_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
    end
    
    % integral initial condition
    y_int = buf_intY(:,max(1,current_pos-1)); 
    % integral
    for i = 1:DynOpt.dim_out
        y_int(i) = trapz(DynOpt.tspan,[y_int(i), y_read(i)]);
    end
end

%%%% 1 OUTPUT %%%%
try 
    if(forward == 1)
        yhat = [x_read, dy, y_int];
    else
        yhat = [x_read, -dy, y_int];
    end
catch
   disp('ARARMAX');
end

%%%% Buffer output %%%%
buf_dy_out = buf_dy;

end