%% signal information policy
function dJ_cond(theta,beta,gamma)

    global DynOpt params
    
    % propagation of Yhat 
    measure_forward = 1;
    buf_temp = DynOpt.buf_dyhat;
    [~, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,buf_temp,DynOpt.intYhat_full_story,params,DynOpt.ActualTimeIndex);
    
    % temp buffers
    Y_temp = DynOpt.Y_full_story;
    dY_temp = DynOpt.dY_full_story;
    Yint_temp = DynOpt.intY_full_story;
    
    Yhat_temp = [DynOpt.Yhat_full_story, Yhat(:,1)];
    dYhat_temp = [DynOpt.dYhat_full_story, Yhat(:,2)];
    Yint_hat_temp = [DynOpt.intYhat_full_story, Yhat(:,3)];
    
    % condition with error
    w = 5;
    mode = 'back';
    pos = size(Y_temp,2);
    DynOpt.e = abs(mean(moving_average_single(Y_temp-Yhat_temp,w,pos,mode)));
    DynOpt.e_dot = abs(mean(moving_average_single(dY_temp-dYhat_temp,w,pos,mode)));
    
    initial_pos = min(w,size(Yint_temp,2)-1);
    temp_e_int = zeros(1,DynOpt.dim_out);
    for i=1:DynOpt.dim_out
        temp_e_int(i) = trapz(DynOpt.Ts,(Yint_temp(i,end-initial_pos:end)-Yint_hat_temp(i,end-initial_pos:end)));
    end
    DynOpt.e_int = abs(mean(temp_e_int));
    
    DynOpt.e_mesh = theta*DynOpt.e + beta*DynOpt.e_dot + gamma*DynOpt.e_int;
    DynOpt.dJ_cond = mean(abs(DynOpt.e_mesh));
    DynOpt.dJ_cond_story(:,end+1) = [DynOpt.e; DynOpt.e_dot; DynOpt.e_int; DynOpt.dJ_cond];
end