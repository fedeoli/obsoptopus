function [buf_dy_out,yhat] = get_measure_v1(x,j,forward,buf_dy)

% Output function to retrieve the estimated output from the estimated state
% NOT CONSIDERING MAGNETOMETERS (only numerical test)
%%%%%%%%%%%%%%%

global DynOpt params 

x_propagate = x;

% check if measure or propagate
if j == 0
    measure_flag = 1;
    time_instant = DynOpt.ActualTimeIndex;
else
    measure_flag = 0;
    time_instant = DynOpt.BackTimeIndex;
end

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously
if (forward==1)
    for i=1:((DynOpt.Nts)*(j-1))    
        time_instant = DynOpt.BackTimeIndex+i;
        x_propagate = model_propagate(time_instant,DynOpt.Ts,x_propagate);  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update buffer
        x_read = zeros(DynOpt.dim_out,1);
        for j=1:length(params.observed_state)
            x_read(j) = x_propagate(params.observed_state(j));
        end
        % derivatives
        dy = zeros(DynOpt.dim_out,1);
        for j=1:DynOpt.dim_out
            [buf_dy(j,:), dy(j)] = IterativePseudoDerivative(DynOpt.Ts,x_read(j),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(j,:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
else
    for i=1:((DynOpt.Nts)*(j-1))
        time_instant = DynOpt.BackTimeIndex-i+1;
        x_propagate = model_propagate(time_instant,DynOpt.Ts,x_propagate);  
    end
end

if j==1 || measure_flag == 1
    %%%% FIRST OBSERVATION
    x_read = zeros(DynOpt.dim_out,1);
    for i=1:length(params.observed_state)
        x_read(i) = x_propagate(params.observed_state(i));
    end 
    % derivatives
    dy = zeros(DynOpt.dim_out,1);
    for j=1:DynOpt.dim_out
        [buf_dy(j,:), dy(j)] = IterativePseudoDerivative(DynOpt.Ts,x_read(j),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(j,:));
    end
end

%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [x_read, dy];
else
    yhat = [x_read, dy];
end

%%%% Buffer output %%%%
buf_dy_out = buf_dy;

end