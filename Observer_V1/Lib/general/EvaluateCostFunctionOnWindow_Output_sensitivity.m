function yhat = EvaluateCostFunctionOnWindow_Output_sensitivity(x,j,forward)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%

global DynOpt params

x_propagate = x(1:DynOpt.StateDim);

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously
% could be speeded up by using (j-1)*Nts:j*Nts or something similar
if (forward==1)
    for i=1:((DynOpt.Nts)*(j-1))
        % actual propagation
        set_input(DynOpt.BackTimeIndex+i);
        x_propagate = PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,1,params);
    end
else
    for i=1:((DynOpt.Nts)*(j-1))
        set_input(DynOpt.BackTimeIndex-i+1);
        x_propagate = PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,-1,params);
    end
end



%%%% FIRST OBSERVATION
x_read(1) = x_propagate(params.observed_state);

% derivative estimate
[DynOpt.buf_dyhat, dx_read] = IterativePseudoDerivative(DynOpt.Ts,x_read(1),DynOpt.c1_derivative,DynOpt.d1_derivative,0,DynOpt.buf_dyhat);

if forward == 1
    yhat = [x_read; dx_read];
else
    yhat = [x_read; -dx_read];
end

yhat = [yhat; zeros(DynOpt.dim_out-length(yhat),1)];

end