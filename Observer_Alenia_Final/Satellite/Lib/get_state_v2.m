%% get state from DynOpt 
function nu = get_state_v2(DynOpt)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION WITH NORM %%%%%%%%%%%%%%%%%%%%%%%
%     % get data
%     e = norm(DynOpt.Yhat_full_story(1:3,end)-DynOpt.Y_full_story(1:3,end));
%     edot = norm(DynOpt.dYhat_full_story(1:3,end)-DynOpt.dY_full_story(1:3,end));
%     
%     % track err
%     e_track = norm(DynOpt.track_err(:,end));
%     edot_track = norm(DynOpt.track_errdot(:,end));
%     
%     % e final
%     theta = 1;
%     e_fin = theta*e + (1-theta)*e_track;
%     theta_dot = 1;
%     edot_fin = theta_dot*edot + (1-theta_dot)*edot_track;
%     
%     % get state
%     nu = [e_fin, edot_fin];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION WITH COMPONENTS %%%%%%%%%%%%%%%%%
%     % get data
%     e = DynOpt.Yhat_full_story(1:3,end)-DynOpt.Y_full_story(1:3,end);
%     edot = DynOpt.dYhat_full_story(1:3,end)-DynOpt.dY_full_story(1:3,end);
%     
%     % track err
%     e_track = DynOpt.track_err(:,end);
%     edot_track = DynOpt.track_errdot(:,end);
%     
%     % e final
%     theta = 1;
%     e_fin = theta*e + (1-theta)*e_track;
%     theta_dot = 1;
%     edot_fin = theta_dot*edot + (1-theta_dot)*edot_track;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION WITH MAGNETOMETERS %%%%%%%%%%%%%%
    % get data
    e_w = DynOpt.Yhat_full_story(1:3,end)-DynOpt.Y_full_story(1:3,end);
    e_mag = DynOpt.Yhat_full_story(4:6,end)-DynOpt.Y_full_story(4:6,end);
    
    % e final
    e_fin = norm(e_w);
    edot_fin = norm(e_mag);
    
    % get state
    nu = [e_fin; edot_fin];
end