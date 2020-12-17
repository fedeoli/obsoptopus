function [buf_dy_out, buf_ddy_out, buf_yint_out, yhat] = EvaluateCostFunctionOnWindow_Output_general_notime_nobuf_params(x,j,forward,buf_dy,buf_ddy,buf_yint)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
% customize code below
%%%%%%%%%%%%%%%

global DynOpt params

x_propagate = x;

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously
% could be speeded up by using (j-1)*Nts:j*Nts or something similar
if (forward==1)
    for i=1:((DynOpt.Nts)*(j-1))       
        set_input(DynOpt.BackTimeIndex+i);
        if strcmp(DynOpt.int_flag,'rk4') 
            tspan = [DynOpt.time(DynOpt.BackTimeIndex+i), DynOpt.time(DynOpt.BackTimeIndex+i)+DynOpt.Ts];
            temp = rk4_V1_1(DynOpt.model, tspan, x_propagate, params);
            x_propagate = temp(:,end);
        else
            x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,1,params);
        end
    end
else
    for i=1:((DynOpt.Nts)*(j-1))
        set_input(DynOpt.BackTimeIndex-i+1);
        if strcmp(DynOpt.int_flag,'rk4') 
            tspan = [DynOpt.time(DynOpt.BackTimeIndex-i+1), DynOpt.time(DynOpt.BackTimeIndex-i+1)-DynOpt.Ts];
            temp = rk4_V1_1(DynOpt.model, tspan, x_propagate, params);
            x_propagate = temp(:,end);
        else
            x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,-1,params);
        end
    end
end



%%%% FIRST OBSERVATION
x_read(1) = x_propagate(params.observed_state)^DynOpt.measure_exp;

% integral action
buf_yint = [buf_yint(1:end-1) x_read(1)];
time_array = 0:DynOpt.Ts:(DynOpt.Ts*(DynOpt.w-1));
y_int(1) = trapz(time_array,buf_yint);

% reset buffer at first iteration
[buf_dy, dy(1)] = IterativePseudoDerivative(DynOpt.Ts,x_read(1),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy);
[buf_ddy, ddy(1)] = IterativePseudoDerivative(DynOpt.Ts,dy(1),DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,buf_ddy);

%%%% 1 OUTPUT %%%%
if(forward == 1)
    yhat = [x_read(1); dy(1); ddy(1); y_int(1)];
else
    yhat = [x_read(1); -dy(1); -ddy(1); -y_int(1)];
end

buf_dy_out = buf_dy;
buf_ddy_out = buf_ddy;
buf_yint_out = buf_yint;

end