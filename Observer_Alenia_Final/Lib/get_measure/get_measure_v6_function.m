%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = get_measure_v6_function(DynOpt,x,j,forward,buf_dy,buf_intY,params,current_pos)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

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
    if DynOpt.bias_enable == 0 || DynOpt.optimise_params == 0
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
    zero_flag = pos(1)>=1;
    % set n_iter
    temp_pos = min(length(pos),j);
    n_iter = sum(buf_dist(1:pos(temp_pos)))-zero_flag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % integral initial condition
    y_int = buf_intY(:,DynOpt.BackIterIndex);
    
    if (forward==1)
        % the max has been added to at least save the current state if
        for i=1:n_iter
            
            % integration
            [xdot,~] = DynOpt.model(x_propagate,params,DynOpt);
            [x_propagate, params] = DynOpt.model_propagate(DynOpt,DynOpt.BackIterIndex+i,DynOpt.Ts,x_propagate,params);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update buffer
            x_read = zeros(9,1);
            for k=1:length(params.observed_state)
                x_read(k) = x_propagate(params.observed_state(k));
            end
                        
            % get magnetometers
            [MagTemp_body1, MagTemp_body2,DynOpt] = Magnetometers_GetMeasures_v3(DynOpt,DynOpt.BackIterIndex+i,x_propagate(offset+1:offset+4)',DynOpt.RPYbetweenMagSensors');
            Mag = [MagTemp_body1; MagTemp_body2];
            
            % output computation
            x_read(4:end) = Mag;
            y_read = x_read;
            
            % derivatives
            dy = zeros(9,1);
            for k=1:9
                [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,y_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
            end  
            % wdot through the model (for param estimation)
%             dy(1:3) = xdot(5:7);
            
            % integral initial condition
            y_int = zeros(9,1);
            temp_pos = max(1,DynOpt.BackIterIndex+i-1);
            if size(buf_intY,2) >= temp_pos
                temp = buf_intY(:,temp_pos); 
                y_int(1:length(temp)) = temp;
            end
            % integral
            for z = 1:9
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
    x_read = zeros(9,1);
    for i=1:length(params.observed_state)
        x_read(i) = x_propagate(params.observed_state(i));
    end 
    
    % get magnetometers
    [MagTemp_body1, MagTemp_body2,DynOpt] = Magnetometers_GetMeasures_v3(DynOpt,current_pos,x_propagate(offset+1:offset+4)',DynOpt.RPYbetweenMagSensors');  
    Mag = [MagTemp_body1; MagTemp_body2];

    % output computation
    x_read(4:end) = Mag;
    y_read = x_read;
    
    % derivatives
    dy = zeros(9,1);
    for k=1:9
        [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,y_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
    end
    
    % integral initial condition
    y_int = zeros(9,1);
    if ~isempty(buf_intY)
        temp = buf_intY(:,max(1,current_pos-1));
        y_int(1:length(temp)) = temp;
    end
    % integral
    for i = 1:9
        y_int(i) = trapz(DynOpt.tspan,[y_int(i), y_read(i)]);
    end
end

%%%% 1 OUTPUT %%%%
yhat = zeros(9,3);
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