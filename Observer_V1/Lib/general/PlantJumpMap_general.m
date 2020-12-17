function Xnew  = PlantJumpMap_general(X,model,k0,forwardpropagation)
% X: state and parameters
% j: index used to retrieve the data (input) on the moving orizon window
% forwardpropagation: 1 if forward propagation, -1 otherwise
% estimate: 1 if used for the estimated plant, 0 othwerwise
% Example: x = A*x+B*u(j); j-1 is the index used to select input value

global DynOpt

x0 = X(1:DynOpt.StateDim);
t0 = DynOpt.time(k0);

if(forwardpropagation == 1) 
    tspan = [t0 t0+DynOpt.Ts];
%     [~, x] = ode45(model,tspan,x0); 
%     x = feval(model,tspan(2),x0);
    x = x0 + model(tspan(1),x0)*DynOpt.Ts;
else
    tspan = [t0 t0-DynOpt.Ts];
%     [~, x] = ode45(model,tspan,x0);
%     x = feval(model,tspan(2),x0);
    x = x0 - model(tspan(1),x0)*DynOpt.Ts;
end



% Xnew = [x(end,:)';X(DynOpt.StateDim+1:end)];
Xnew = [x; X(DynOpt.StateDim+1:end)];