%% signal information policy
function dJ_cond(theta)

    global DynOpt params
    
    % propagation of Yhat 
    measure_forward = 1;
    buf_temp = DynOpt.buf_dyhat;
    
    % get measure
    [~, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,buf_temp,DynOpt.intYhat_full_story,params);
    
    % temp buffers
    Y_temp = DynOpt.Y_full_story;
    dY_temp = DynOpt.dY_full_story;
    
    Yhat_temp = [DynOpt.Yhat_full_story, Yhat(:,1)];
    dYhat_temp = [DynOpt.dYhat_full_story, Yhat(:,2)];
    
    % condition with error
    w = 5;
    mode = 'back';
    pos = size(Y_temp,2);
    DynOpt.e = moving_average_single(Y_temp-Yhat_temp,w,pos,mode);
    DynOpt.e_dot = moving_average_single(dY_temp-dYhat_temp,w,pos,mode);
    DynOpt.e_mesh = theta*DynOpt.e + (1-theta)*DynOpt.e_dot;
    DynOpt.dJ_cond = abs(DynOpt.e_mesh);
    DynOpt.dJ_cond_story(:,end+1) = [DynOpt.e; DynOpt.e_dot; DynOpt.e_mesh];

end