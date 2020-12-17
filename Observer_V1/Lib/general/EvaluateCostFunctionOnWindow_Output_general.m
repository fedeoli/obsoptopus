function yhat = EvaluateCostFunctionOnWindow_Output_general(x,j,forward)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%
global DynOpt

x_propagate = x;

% start from x0 --> propagate up to x at the jth window
% could be speeded up by using (j-1)*Nts:j*Nts or something similar
for i=1:j*DynOpt.Nts
    if (forward==1)
        x_propagate = PlantJumpMap_general(x_propagate,DynOpt.model,DynOpt.BackTimeIndex+i-1,1);
    else
        x_propagate = PlantJumpMap_general(x_propagate,DynOpt.model,DynOpt.BackTimeIndex-i+1,-1);
    end
end

%%%% FIRST OBSERVATION
x_read(1) = x_propagate(3);

% reset buffer at first iteration
if DynOpt.ActualTimeIndex == 1
    dy(1) = IterativePseudoDerivative_1hat(DynOpt.Ts,x_read(1),DynOpt.c1_derivative,DynOpt.d1_derivative,0,1);
    ddy(1) = IterativePseudoDerivative_2hat(DynOpt.Ts,dy(1),DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,1);
else
    dy(1) = IterativePseudoDerivative_1hat(DynOpt.Ts,x_read(1),DynOpt.c1_derivative,DynOpt.d1_derivative,0,0);
    ddy(1) = IterativePseudoDerivative_2hat(DynOpt.Ts,dy(1),DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,0);    
end

%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [x_read(1); dy(1); ddy(1)];
else
    yhat = [x_read(1); -dy(1); -ddy(1)];
end

end