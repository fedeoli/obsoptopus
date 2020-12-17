function [buf_dy_out, yhat] = EvaluateCostFunctionOnWindow_Output_general_data(j,forward,buf_dy)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%

global DynOpt

%%%% FIRST OBSERVATION
x_read(1) = DynOpt.state(1,j);

% reset buffer at first iteration
[buf_dy, dy] = IterativePseudoDerivative(DynOpt.Ts,x_read(1),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy);


%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [x_read(1), dy];
else
    yhat = [x_read(1), -dy];
end

buf_dy_out = buf_dy;
end