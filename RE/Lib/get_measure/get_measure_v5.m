%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = get_measure_v5(x,j,forward,buf_dy,buf_intY,params,current_pos)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt

x_propagate = x;

offset = DynOpt.integration_pos*6;
if DynOpt.integration_pos
    position = x(1:3)';
else
    position = params.SatellitesCoordinates(1:3)';
end

% check if measure or propagate
if j == 0
    measure_flag = 1;
else
    measure_flag = 0;
    % store bias
    if DynOpt.bias_enable == 0
        b = zeros(DynOpt.nbias,1);
    else
        b = x_propagate(offset+8:offset+8+DynOpt.nbias-1);
    end
    % generate bias
    bias = zeros(DynOpt.dim_out,1);
    for z = 1:length(b)
        bias(z) = bias(z)+b(z);
    end
end

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously. This is done only if measure_flag
% == 0 namely only during optimisation
if measure_flag == 0
    % get adaptive buffer distance. The zero is added to make the vector of the
    % same length of the buffer as the diff command computes ((i+1)th-ith)
    buf_dist = [max(DynOpt.Y_space(1)-DynOpt.BackTimeIndex,0), diff(DynOpt.Y_space)];
    % get first nonzero element
    [~,pos] = find(buf_dist);
    % check if exist a zero
    zero_flag = pos(1)>1;
    % set n_iter
    n_iter = sum(buf_dist(1:pos(j)))-zero_flag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % integral initial condition
    y_int = buf_intY(:,DynOpt.BackTimeIndex);
    
    if (forward==1)
        % the max has been added to at least save the current state if
        for i=1:n_iter
            
            % integration
            [x_propagate, params] = DynOpt.model_propagate(DynOpt.BackTimeIndex+i,DynOpt.Ts,x_propagate,params);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update buffer
            x_read = zeros(DynOpt.dim_out,1);
            for k=1:length(params.observed_state)
                x_read(k) = x_propagate(params.observed_state(k));
            end
                        
            % get magnetometers
            if DynOpt.nMagneto > 0
                [MagTemp_body1, MagTemp_body2,DynOpt] = Magnetometers_GetMeasures_v3(DynOpt,DynOpt.BackTimeIndex+i,x_propagate(offset+1:offset+4)',DynOpt.RPYbetweenMagSensors');
            end
            if DynOpt.nMagneto == 1
                Mag = MagTemp_body1;
            elseif DynOpt.nMagneto == 2
                Mag = [MagTemp_body1; MagTemp_body2];
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
            y_int = buf_intY(:,max(1,DynOpt.BackTimeIndex+i-1)); 
            % integral
            for z = 1:DynOpt.dim_out
                y_int(z) = trapz(DynOpt.tspan,[y_int(z), y_read(z)]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
              
        % edit initial state
        if DynOpt.bias_enable == 1
            x_read = x_read+bias;
        end
    else
        %%%%% TO BE DONE %%%%
    end
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
        [MagTemp_body1, MagTemp_body2,DynOpt] = Magnetometers_GetMeasures_v3(DynOpt,current_pos,x_propagate(offset+1:offset+4)',DynOpt.RPYbetweenMagSensors');
    end
    
    if DynOpt.nMagneto == 1
        Mag = MagTemp_body1;
    elseif DynOpt.nMagneto == 2
        Mag = [MagTemp_body1; MagTemp_body2];
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
    if ~isempty(buf_intY)
        y_int = buf_intY(:,max(1,current_pos-1)); 
    else
        y_int = zeros(DynOpt.dim_out,1);
    end
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