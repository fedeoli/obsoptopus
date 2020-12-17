%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = EvaluateCostFunctionOnWindow_Output_v1(x,j,forward,buf_dy)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
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
        set_input(DynOpt.BackTimeIndex+i);
        x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,1,params);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update buffer
        x_read = zeros(DynOpt.dim_out,1);
        for k=1:length(params.observed_state)
            x_read(k) = x_propagate(params.observed_state(k));
        end
        % derivatives
        dy = zeros(DynOpt.dim_out,1);
        for k=1:DynOpt.dim_out
            [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,x_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
else
    for i=1:((DynOpt.Nts)*(j-1))
        set_input(DynOpt.BackTimeIndex-i+1);
        x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,-1,params);
    end
end

%%%% FIRST OBSERVATION
if j==1 || measure_flag == 1
    %%%% FIRST OBSERVATION
    x_read = zeros(DynOpt.dim_out,1);
    for i=1:length(params.observed_state)
        x_read(i) = x_propagate(params.observed_state(i));
    end 
    % derivatives
    dy = zeros(DynOpt.dim_out,1);
    for k=1:DynOpt.dim_out
        [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,x_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
    end
end

%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [x_read, dy];
else
    yhat = [x_read, -dy];
end

%%%% Buffer output %%%%
buf_dy_out = buf_dy;

end