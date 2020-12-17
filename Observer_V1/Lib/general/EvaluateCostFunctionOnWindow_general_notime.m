function Jtot = EvaluateCostFunctionOnWindow_general_notime(x)

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

%%%%% WORKSPACE %%%%%
% DynOpt.ActualTimeIndex-(DynOpt.w-j)*DynOpt.Nts

for j=1:DynOpt.w
    
    %evaluate the weighted in time cost function at this iteration time
    if(DynOpt.ForwardOptimization ~= 1) %backward
        
        Yhat = EvaluateCostFunctionOnWindow_Output_general_notime(X,j-1,-1);
        
        Jtot = Jtot + DynOpt.Weight(end-j+1)*(DynOpt.Y(1,end-j+1)-Yhat(1))^2 + ... 
                           DynOpt.dWeight(end-j+1)*(DynOpt.Y(2,end-j+1)-Yhat(2))^2+...
                           DynOpt.ddWeight(end-j+1)*(DynOpt.Y(3,end-j+1)-Yhat(3))^2;
                       
%%% 2 OUTPUTS COST FUNCTION                       
%                            DynOpt.Weight(end-j+1)*(DynOpt.Y(4,end-j+1)-Yhat(4))^2 + ... 
%                            DynOpt.dWeight(end-j+1)*(DynOpt.Y(5,end-j+1)-Yhat(5))^2+...
%                            DynOpt.ddWeight(end-j+1)*(DynOpt.Y(6,end-j+1)-Yhat(6))^2;                       
        
    else %forward    k-(DynOpt.w)*DynOpt.Nts+1
        
        Yhat = EvaluateCostFunctionOnWindow_Output_general_notime(X,j-1,1);
        
        Jtot = Jtot + DynOpt.Weight(j)*(DynOpt.Y(1,j)-Yhat(1))^2 + ... 
                           DynOpt.dWeight(j)*(DynOpt.Y(2,j)-Yhat(2))^2+...
                           DynOpt.ddWeight(j)*(DynOpt.Y(3,j)-Yhat(3))^2;
                       
%%% 2 OUTPUTS COST FUNCTION                       
%                            DynOpt.Weight(end-j+1)*(DynOpt.Y(4,end-j+1)-Yhat(4))^2 + ... 
%                            DynOpt.dWeight(end-j+1)*(DynOpt.Y(5,end-j+1)-Yhat(5))^2+...
%                            DynOpt.ddWeight(end-j+1)*(DynOpt.Y(6,end-j+1)-Yhat(6))^2; 
                       
    end
    
end

%Jtot
    