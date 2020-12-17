function Jtot = EvaluateCostFunctionOnWindow_general_notime_nobuf(x)

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
        
        [DynOpt.buf_dyhat,DynOpt.buf_ddyhat,Yhat] = EvaluateCostFunctionOnWindow_Output_general_notime_nobuf(X,j,-1,DynOpt.buf_dyhat,DynOpt.buf_ddyhat);
        
        % moving average - filter data to increase derivatives estimation
        if DynOpt.filter_flag_opt == 1 && DynOpt.ActualTimeIndex > DynOpt.filter_window 
            Y_buf = [DynOpt.Y_full_story(:,1:DynOpt.ActualTimeIndex-1), Yhat];
            pos = DynOpt.ActualTimeIndex;
            temp = moving_average(Y_buf,DynOpt.filter_window_opt,pos,DynOpt.filter_mode);
            Y_filter_hat = temp(:,pos);
        else
            Y_filter_hat = Yhat;
        end
        
        Jtot = Jtot + DynOpt.Weight(end-j)*(DynOpt.Y(1,end-j)-Y_filter_hat(1))^2 + ... 
                           DynOpt.dWeight(end-j)*(DynOpt.Y(2,end-j+1)-Y_filter_hat(2))^2 + ...
                           DynOpt.ddWeight(end-j)*(DynOpt.Y(3,end-j+1)-Y_filter_hat(3))^2;                                        
        
    else %forward    
        
        for i=(j-1):-1:1
            DynOpt.buf_dyhat(end-i+1) = DynOpt.buf_dy(end-DynOpt.w + (j-i));
            DynOpt.buf_ddyhat(end-i+1) = DynOpt.buf_ddy(end-DynOpt.w + (j-i));
        end
        
        [DynOpt.buf_dyhat,DynOpt.buf_ddyhat,Yhat] = EvaluateCostFunctionOnWindow_Output_general_notime_nobuf(X,j,1,DynOpt.buf_dyhat,DynOpt.buf_ddyhat);
        
        % moving average - filter data to increase derivatives estimation
        if DynOpt.filter_flag_opt == 1 && DynOpt.ActualTimeIndex > DynOpt.filter_window_opt 
            Y_buf = [DynOpt.Y_full_story(:,1:DynOpt.ActualTimeIndex-1), Yhat];
            pos = length(Y_buf);
            temp = moving_average(Y_buf,DynOpt.filter_window_opt,pos,DynOpt.filter_mode);
            Y_filter_hat = temp(:,pos);
        else
            Y_filter_hat = Yhat;
        end
        
        Jtot = Jtot + DynOpt.Weight(j)*(DynOpt.Y(1,j)-Y_filter_hat(1))^2 + ... 
                           DynOpt.dWeight(j)*(DynOpt.Y(2,j)-Y_filter_hat(2))^2 + ...
                           DynOpt.ddWeight(j)*(DynOpt.Y(3,j)-Y_filter_hat(3))^2;
                       
    end
    
end
    