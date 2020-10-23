function yhat = EvaluateCostFunctionOnWindow_EstimatedOutput_v1(x)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%

%x = [int(z), z, dz, I_H, beta, gamma_0, gamma_1]'
%u = [I_href, elong]
yhat = x(1);%x(2);x(3),x(4)];
