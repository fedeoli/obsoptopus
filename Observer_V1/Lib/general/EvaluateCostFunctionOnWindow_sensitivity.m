function Jtot = EvaluateCostFunctionOnWindow_sensitivity(x)

global DynOpt

%x: the estimated optimal state and parameters vector provided by fmin...

%DynOpt.Ts : Sampling Time
%DynOpt.Nts : Inter-sampling Time
%DynOpt.X : actual state (t_k), a column vector, first the state, then stacked the parameters
%DynOpt.U : input vector, last value at time t_k, a matrix with w column vectors with the number of rows as the number of "instantaneous" input
%DynOpt.Y : output vector, last value at time t_k, a matrix with w column vectors with the number of rows as the number of "instantaneous" outputs
%DynOpt.w : window size
%DynOpt.StateDim : number of state variables (subject to dynamics, not parameters)
%DynOpt.Weight : window weight size, last on the more recent t_k

% cost function init
Jtot = 0;

%optimization vector  
X = x; 

% dyhat
windowsamples = floor(DynOpt.WindowSamples/(DynOpt.Nts+1));
checkpoint = max(0,floor(DynOpt.ActualTimeIndex/(DynOpt.Nts+1))-windowsamples);

if checkpoint >= DynOpt.d1_derivative
    checkpoint = DynOpt.d1_derivative;
end

if checkpoint > 0
   nzeros = DynOpt.d1_derivative-checkpoint;
   range = DynOpt.BackTimeIndex-DynOpt.Nts*checkpoint:DynOpt.Nts:DynOpt.BackTimeIndex-DynOpt.Nts;
   DynOpt.buf_dyhat = [zeros(1,nzeros), DynOpt.OptXstory(DynOpt.params.observed_state,range)];
end

for j=1:DynOpt.w
    
    %evaluate the weighted in time cost function at this iteration time
    if(DynOpt.ForwardOptimization ~= 1) %backward
        
        Yhat = EvaluateCostFunctionOnWindow_Output_sensitivity(X,j,-1);
        
        J1 = (norm(DynOpt.Y(:,end-j+1)-Yhat))^2;
        
        Jtot = Jtot + J1;
       
    else %forward 
        
        Yhat = EvaluateCostFunctionOnWindow_Output_sensitivity(X,j,1);
        
        J1 = (norm(DynOpt.Y(1:DynOpt.y_end,j)-Yhat(1:DynOpt.y_end)))^2;
        
        Jtot = Jtot + J1;
                       
    end
    
end

end
    