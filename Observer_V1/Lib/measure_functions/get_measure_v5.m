%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = get_measure_v5(x,j,forward,buf_dy,buf_intY,params)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt

x_propagate = x;

% check if measure or propagate
if j == 0
    measure_flag = 1;
else
    measure_flag = 0;
end

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously. This is done only if measure_flag
% == 0 namely only during optimisation
if measure_flag == 0
    %%%%%%%%%%%%%%% DEFINE NUMBER OF INTEGRATIONS %%%%%%%%%%%%%%%
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
    
    % integration
    if (forward==1)
        % the max has been added to at least save the current state if
        % diff(j) == 0 (no propagation at all)
        for i=1:n_iter
            
            set_input(DynOpt.BackTimeIndex+i);
            temp_state = rk4_V1_1(DynOpt.model, DynOpt.tspan, x_propagate, params);   
            x_propagate =  temp_state(:,end);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update buffer
            x_read = zeros(DynOpt.dim_out,1);
            for k=1:length(params.observed_state)
                x_read(k) = x_propagate(params.observed_state(k));
            end
            
            % compute output
            y_read = x_read(1)+abs(x_read(2));
            
            % derivatives
            dy = zeros(DynOpt.dim_out,1);
            for k=1:DynOpt.dim_out
                [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,y_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
            end
            
            % integral
            y_int = trapz(DynOpt.tspan,[y_int, y_read]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
    else
        %%%%% TO BE DONE %%%%
    end
    
    
    
end

%%%% FIRST OBSERVATION
if measure_flag == 1
    % state
    x_read = zeros(DynOpt.dim_out,1);
    for i=1:length(params.observed_state)
        x_read(i) = x_propagate(params.observed_state(i));
    end 
    
    % compute output
    y_read = x_read(1)+abs(x_read(2));
            
    % derivatives
    dy = zeros(DynOpt.dim_out,1);
    for k=1:DynOpt.dim_out
        [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,y_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
    end
    
    % integral initial condition
    if DynOpt.synthetic_flag == 0
       y_int = buf_intY(:,max(1,DynOpt.ActualTimeIndex-1)); 
    else
       y_int = 0; 
    end
    % integral
    y_int = trapz(DynOpt.tspan,[y_int, y_read]);
end

%%%% 1 OUTPUT %%%%
try 
    if(forward == 1)
        yhat = [y_read, dy, y_int];
    else
        yhat = [y_read, -dy, y_int];
    end
catch
   disp('ARARMAX');
end

%%%% Buffer output %%%%
buf_dy_out = buf_dy;

end