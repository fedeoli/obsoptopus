function Xfinal = EvaluateCostFunctionOnWindow_StateBackInTimePropagation_v1(X,i,forward)
%% Integration function of the dynamics forward (+1) or backward (-1) in time between intersamples
%X : state and parameters
%i : intersampled counter
%forward: (+1) integration forward in time, (-1) backward in time.
global DynOpt
x = X(1:DynOpt.StateDim);
gamma0_hat = X(5);
gamma1_hat = X(6);
for k = (i-1)*DynOpt.Nts+1:(i)*DynOpt.Nts,
    A = [0,1,0,0; 0,0,1,0; 0,0,0,gamma0_hat+gamma1_hat*(DynOpt.U(2,end-k)-1);-DynOpt.Ki*DynOpt.beta,-DynOpt.Kp*DynOpt.beta,-DynOpt.Kd*DynOpt.beta,-DynOpt.beta];
    x = x + forward*DynOpt.Ts*(A*x);
end
Xfinal = x;

