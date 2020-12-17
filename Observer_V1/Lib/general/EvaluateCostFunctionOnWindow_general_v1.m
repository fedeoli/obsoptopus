%% cost function
function Jtot = EvaluateCostFunctionOnWindow_general_v1(x)

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

% params update
params_update(X);

% set the derivative buffer as before the optimisation process (multiple f computation)
if DynOpt.BackTimeIndex >= DynOpt.d1_derivative
    DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.d1_derivative-1):DynOpt.BackTimeIndex);
else
    init_pos = DynOpt.d1_derivative-DynOpt.BackTimeIndex;
    DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.BackTimeIndex-1):DynOpt.BackTimeIndex)];
end

for j=1:DynOpt.w
    
    %evaluate the weighted in time cost function at this iteration time
    if(DynOpt.ForwardOptimization ~= 1) %backward        
        %%%%%% TO BE DONE %%%%%%
    else  
        
        [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(X,j,1,DynOpt.buf_dyhat_temp);
                
        J_meas = zeros(1,size(Yhat,1));
        J_der = zeros(1,size(Yhat,1));
        
        for i=1:length(J_meas)
            diff = (DynOpt.Y(i,j)-Yhat(i,1));
            J_meas(i) = DynOpt.scale_factor(1,i)*(diff)^2;
        end
        
        % not considering magnetometers 
        n_der = DynOpt.dim_out;
        for i=1:n_der
            diff = (DynOpt.dY(i,j)-Yhat(i,2));
            J_der(i) = DynOpt.scale_factor(2,i)*(diff)^2;
        end
        
        Jtot = Jtot + sum(J_meas) + sum(J_der);     
        
        % save J_meas and J_der
        DynOpt.J_meas(end+1) = sum(J_meas);
        DynOpt.J_der(end+1) = sum(J_der);
                       
    end
    
end
    DynOpt.Yhat_temp = Yhat;
end
    