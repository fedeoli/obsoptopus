function yhat = EvaluateCostFunctionOnWindow_Output_v1(x,k,forward)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%
global DynOpt

%x = [int(z), z, dz, I_H, beta, gamma_0, gamma_1]'
%u = [I_href, elong]
dy = IterativePseudoDerivative_1hat(DynOpt.Ts,x(2),DynOpt.c1_derivative,DynOpt.d1_derivative,0,0);
ddy = IterativePseudoDerivative_2hat(DynOpt.Ts,dy,DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,0);
if(forward == 1)
    yhat = [x(2);dy;ddy];
else
    yhat = [x(2);-dy;-ddy];
end