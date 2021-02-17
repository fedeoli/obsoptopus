%% noise model - gyroscope in this case
function Y_noise = noise_model_v1(Y_true,DynOpt,params)

    b = DynOpt.OptXstoryTRUE(8,DynOpt.ActualTimeIndex);
    eps_noise = [DynOpt.measure_noise.*ones(length(params.observed_state),1); zeros(DynOpt.dim_out-length(params.observed_state),1)];
    
    % copy true on derivative and integral
    Y_noise = Y_true;
    
    % edit the data measurement
    Y_noise(:,1) = Y_true(:,1) + b + eps_noise;

end