function Xnew  = PlantJumpMap_general_notime_params_hybrid(X1,X2,model,forwardpropagation,params)
% X: state and parameters
% j: index used to retrieve the data (input) on the moving orizon window
% forwardpropagation: 1 if forward propagation, -1 otherwise
% estimate: 1 if used for the estimated plant, 0 othwerwise
% Example: x = A*x+B*u(j); j-1 is the index used to select input value

global DynOpt

x01 = X1(1:DynOpt.StateDim+DynOpt.nparams);
x02 = X2;

tjunk = 0;

if(forwardpropagation == 1) 
    x = x02 + model(tjunk,x01,x02,params)*DynOpt.Ts;
else
    x = x02 - model(tjunk,x01,x02,params)*DynOpt.Ts;
end



% Xnew = [x(end,:)';X(DynOpt.StateDim+1:end)];
Xnew = x;