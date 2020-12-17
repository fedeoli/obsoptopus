function [buf_dy_out, buf_ddy_out, yhat] = EvaluateCostFunctionOnWindow_Output_general_notime_nobuf_v2(x,j,forward,buf_dy,buf_ddy)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%

global DynOpt params

x_propagate = x;

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously
% could be speeded up by using (j-1)*Nts:j*Nts or something similar
for i=1:j*DynOpt.Nts
    if (forward==1)
        x_propagate = PlantJumpMap_general_notime(x_propagate,DynOpt.model,1);
    else
        x_propagate = PlantJumpMap_general_notime(x_propagate,DynOpt.model,-1);
    end
end

%%%% FIRST OBSERVATION
x_read(1) = x_propagate(params.observed_state);

% substitute if zero vals detected
pos = min(DynOpt.d1_derivative,DynOpt.ActualTimeIndex);
for i=1:pos
   if buf_dy(i) == 0
      buf_dy(i) = DynOpt.Y_full_story(2,DynOpt.ActualTimeIndex-pos+i); 
   end
end

pos = min(DynOpt.d1_dderivative,DynOpt.ActualTimeIndex);
for i=1:pos
   if buf_ddy(i) == 0
      buf_ddy(i) = DynOpt.Y_full_story(3,DynOpt.ActualTimeIndex-pos+i); 
   end 
end

% reset buffer at first iteration
[buf_dy, dy(1)] = IterativePseudoDerivative(DynOpt.Ts,x_read(1),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy);
[buf_ddy, ddy(1)] = IterativePseudoDerivative(DynOpt.Ts,dy(1),DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,buf_ddy);

%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [x_read(1); dy(1); ddy(1)];
else
    yhat = [x_read(1); -dy(1); -ddy(1)];
end

buf_dy_out = buf_dy;
buf_ddy_out = buf_ddy;

end