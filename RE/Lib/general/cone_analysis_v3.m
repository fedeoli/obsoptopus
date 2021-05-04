%% define and plot cone analysis on the estimated model
%%% global init %%%
global DynOpt
% free to upsate the input
DynOpt.recollect_input = 1;

% set synthetic data flag
DynOpt.synthetic_int = 1;

% scale up
scale = 5;

% clear 
clear x_start_low_story
clear x_start_high_story
clear DynOpt.coneState_low
clear DynOpt.coneState_high
clear DynOpt.cone_quat_low
clear DynOpt.cone_quat_high

for pos = 2:DynOpt.Niter

    %%%%% init section %%%%%
    time = max(1,DynOpt.time(pos));
    start_pos = max(1,pos-5);
    time_array = DynOpt.time(start_pos:pos);
    offset = DynOpt.integration_pos*6;
    

    % temp buffers
    Y_temp = DynOpt.Y_full_story(:,start_pos:pos);
    Yhat_temp = DynOpt.Yhat_full_story(:,start_pos:pos);

    %%% sigma - angular velocity %%%
    int_arg = 1/time*(Y_temp(1:3,:)-Yhat_temp(1:3,:)).^2;
    sigma_w = zeros(3,1);
    for i=1:3
        sigma_w(i) = scale*trapz(time_array,int_arg(i,:));
    end

    %%% sigma - quaternion %%%
    int_arg = 1/time*(Y_temp(4:end,:)-Yhat_temp(4:end,:)).^2;
    sigma_eul = zeros(3,1);
    for i=1:3
        sigma_eul(i) = scale*trapz(time_array,int_arg(i,:));
    end

    % store sigma
    DynOpt.sigma_w(:,pos) = sigma_w;
    DynOpt.sigma_eul(:,pos) = sigma_eul;   

    % init state - low
    x_start_low = DynOpt.OptXstory_runtime(:,pos);
    quat_start = x_start_low(offset+1:offset+4);
    eul_start = quat2eul(transpose(quat_start));
    x_start_low(offset+1:offset+4) = eul2quat(eul_start-transpose(sigma_eul));
    x_start_low(offset+5:offset+7) = x_start_low(offset+5:offset+7)-sigma_w;
    x_start_low_story(:,pos) = x_start_low;
    

    % init state - high
    x_start_high = DynOpt.OptXstory_runtime(:,pos);
    quat_start = x_start_high(offset+1:offset+4);
    eul_start = quat2eul(transpose(quat_start));
    x_start_high(offset+1:offset+4) = eul2quat(eul_start+transpose(sigma_eul));
    x_start_high(offset+5:offset+7) = x_start_high(offset+5:offset+7)+sigma_w;
    x_start_high_story(:,pos) = x_start_high;


    % integration - low
%     [x_low_story(:,pos-1), params] = DynOpt.model_propagate(i,DynOpt.Ts,x_start_low,params);
%     % integration - high
%     [x_high_story(:,pos-1), params] = DynOpt.model_propagate(i,DynOpt.Ts,x_start_high,params);

end

% state storage
DynOpt.coneState_low = x_start_low_story;
DynOpt.coneState_high = x_start_high_story;


% store quaternions in euler angles
DynOpt.cone_quat_low = DynOpt.wrap(quat2eul(DynOpt.coneState_low(offset+1:offset+4,:)'))';
DynOpt.cone_quat_high = DynOpt.wrap(quat2eul(DynOpt.coneState_high(offset+1:offset+4,:)'))';

% set synthetic data flag
DynOpt.synthetic_int = 0;

% clear temp vars
keep DynOpt