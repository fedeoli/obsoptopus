%% noise model - gyroscope in this case
function Y_noise = noise_model_v2(Y_true,DynOpt,params)

    % bias
    if DynOpt.bias_dyn == 1
        b = DynOpt.param_story;
    else
        b = 0;
    end
    
    % measurement noise
    DynOpt.measure_noise = DynOpt.measure_amp*randn(length(params.observed_state),DynOpt.Niter);
    eps_noise = [DynOpt.measure_noise.*ones(length(params.observed_state),1); zeros(DynOpt.dim_out-length(params.observed_state),1)];
    
    % edit the data measurement
    Y_noise = Y_true + b + eps_noise;

end