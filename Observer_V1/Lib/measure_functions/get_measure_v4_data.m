function [buf_dy_out, yhat] = get_measure_v4_data(j,forward,buf_dy)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%

global DynOpt 

% compute output
y_read = DynOpt.stateStory(1,j);

% reset buffer at first iteration
[buf_dy, dy] = IterativePseudoDerivative(DynOpt.Ts,y_read,DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy);


%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [y_read, dy];
else
    yhat = [y_read, -dy];
end

buf_dy_out = buf_dy;
end