%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = EvaluateCostFunctionOnWindow_Output_v2(x,j,forward,buf_dy)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt params

x_propagate = x;

% check if measure or propagate
if j == 0
    measure_flag = 1;
else
    measure_flag = 0;
end

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously
if (forward==1)
    for i=1:((DynOpt.Nts)*(j-1))       
        set_input(DynOpt.BackTimeIndex+i);
        temp_state = rk4_V1_1(DynOpt.model, DynOpt.tspan, x_propagate, params);   
        x_propagate =  temp_state(:,end);
        
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
    %%%%% TO BE DONE %%%%
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