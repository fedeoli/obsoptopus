%% signal information policy
function dJ_cond_sensitivity_v1(theta)

    global DynOpt
    
    % propagation of Yhat 
    measure_forward = 1;
    buf_temp = DynOpt.buf_dyhat;
    buf2_temp = DynOpt.buf_ddyhat;
    [~, ~, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,buf_temp,buf2_temp);
    
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
    DynOpt.dJ_cond = abs(DynOpt.e_mesh) <= DynOpt.dJcond_thresh;
    DynOpt.dJ_cond_story(:,end+1) = [DynOpt.e; DynOpt.e_dot; DynOpt.e_mesh];
    
    % condition with J function
%     DynOpt.e_mesh = theta*DynOpt.Jstory(end) + (1-theta)*DynOpt.Jdot_story(end);
%     DynOpt.dJ_cond = (abs(DynOpt.e_mesh) <= DynOpt.dJ_thresh) && (DynOpt.ActualTimeIndex > DynOpt.WindowSamples);
%     DynOpt.dJ_cond_story(:,end+1) = [DynOpt.Jstory(end); DynOpt.Jdot_story(end); DynOpt.e_mesh];
    
    % condition with simple output
%         dy_small = (abs(DynOpt.dY_full_story(:,end)) <= DynOpt.dy_thresh);

end