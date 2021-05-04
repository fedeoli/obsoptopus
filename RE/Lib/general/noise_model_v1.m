%% noise model - gyroscope in this case
function Y_noise = noise_model_v1(Y_true,DynOpt,params)
    
    % offset
    offset = DynOpt.integration_pos*6;

    % bias
    if DynOpt.bias_enable == 1
        b = DynOpt.OptXstoryTRUE(offset+8,DynOpt.ActualTimeIndex);
    else
        b = 0;
    end
    
    % measurement noise
    DynOpt.measure_noise = DynOpt.eps_noise_story(:,DynOpt.ActualTimeIndex);
%     eps_noise = [DynOpt.measure_noise.*ones(length(params.observed_state),1); zeros(DynOpt.dim_out-length(params.observed_state),1)];
    eps_noise = DynOpt.measure_noise;
    
    % copy true on derivative and integral
    Y_noise = Y_true;
    
    % edit the data measurement
    gyrobias = b*ones(3,1);
    magbias = params.MagnetoBias(1:3*DynOpt.nMagneto);
    bias = [gyrobias;magbias];
    Y_noise(:,1) = Y_true(:,1) + bias + eps_noise;

end