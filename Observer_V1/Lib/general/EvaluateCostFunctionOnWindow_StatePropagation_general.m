function Xfinal = EvaluateCostFunctionOnWindow_StatePropagation_general(X,j)
%% Integration function of the dynamics forward (+1) or backward (-1) in time between intersamples
%X : state and parameters
%j : intersampled counter of the window 

global DynOpt

if(DynOpt.ForwardOptimization ~= 1)
    %integrating backward between intersamples
    for k = (j-1)*DynOpt.Nts+1:(j)*DynOpt.Nts %use it as: x(end-k-1) = x(end-k) + Ts*(A*x(end-k)+Bu(end-k)).
        X = PlantJumpMap_general(X,DynOpt.model,k,-1);
    end
else
    %forward integration
    for k = (j-1)*DynOpt.Nts+1:(j)*DynOpt.Nts %use it as: x(k) = x(k-1) + Ts*(A*x(k-1)+Bu(k-1)).
        X = PlantJumpMap_general(X,DynOpt.model,k,1);
    end
end
Xfinal = X;

 

