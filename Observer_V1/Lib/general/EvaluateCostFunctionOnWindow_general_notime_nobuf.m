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
        len_story = length(DynOpt.Y_full_story);
        if DynOpt.filter_flag_opt == 1 && (len_story > DynOpt.filter_window_opt)
            Y_buf = [DynOpt.Y_full_story(:,1:len_story-1), Yhat];
            pos = length(Y_buf);
            temp = moving_average(Y_buf,DynOpt.filter_window_opt,pos,DynOpt.filter_mode);
            Y_filter_hat = temp(:,pos);
        else
            Y_filter_hat = Yhat;
        end
        
        Jtot = Jtot + DynOpt.Weight(end-j+1)*(DynOpt.Y(1,end-j+1)-Y_filter_hat(1))^2 + ... 
                           DynOpt.dWeight(end-j+1)*(DynOpt.Y(2,end-j+1)-Y_filter_hat(2))^2 + ...
                           DynOpt.ddWeight(end-j+1)*(DynOpt.Y(3,end-j+1)-Y_filter_hat(3))^2;                                        
       
    else %forward 
        
        [DynOpt.buf_dyhat,DynOpt.buf_ddyhat,Yhat] = EvaluateCostFunctionOnWindow_Output_general_notime_nobuf(X,j,1,DynOpt.buf_dyhat,DynOpt.buf_ddyhat);
        
        % moving average - filter data to increase derivatives estimation
        len_story = length(DynOpt.Y_full_story);
        if DynOpt.filter_flag_opt == 1 && (len_story > DynOpt.filter_window_opt)
            Y_buf = [DynOpt.Y_full_story(:,1:len_story-1), Yhat];
            pos = length(Y_buf);
            temp = moving_average(Y_buf,DynOpt.filter_window_opt,pos,DynOpt.filter_mode);
            Y_filter_hat = temp(:,pos);
        else
            Y_filter_hat = Yhat;
        end
        
        J1 = (1e0*(DynOpt.Y(1,j)-Y_filter_hat(1)))^2;
        J2 = (1e-3*(DynOpt.Y(2,j)-Y_filter_hat(2)))^2;
        J3 = (DynOpt.Y(3,j)-Y_filter_hat(3))^2;
        Jtot = Jtot + DynOpt.Weight(j)*J1 + ... 
                           DynOpt.dWeight(j)*J2 + ...
                           DynOpt.ddWeight(j)*J3;
                       
    end
    
end
    