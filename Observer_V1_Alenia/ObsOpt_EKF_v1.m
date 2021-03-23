%% ObsOpt implementation
function [DynOpt_out, params_out] = ObsOpt_EKF_v1

global DynOpt params

%%% init EKF %%%
disp('INIT EKF')
[DynOpt,params] = SymAnalysis_RL_v3;
DynOpt.P = 1e-3*eye(7);
DynOpt.Q = 1*1e-6*eye(7);
DynOpt.R = 1*1e-3*eye(DynOpt.dim_out);
DynOpt.recollect_input = 1;

disp('Processing data with the optimization-based observer...')
run_time = tic;
for k=1:length(DynOpt.time)

    % update actual index
    DynOpt.ActualTimeIndex = k;

    % reference state - used for noise
    DynOpt.Xtrue = [DynOpt.state(:,k);DynOpt.param_story(:,k)];

    %forward propagation of the previous estimate
    if(k>1)
        % INTEGRATION OF BOTH POSITION AND ATTITUDE - STACK 
        % Control allocation inside "params" structure          
        [DynOpt.X, params] = DynOpt.model_propagate(DynOpt.ActualTimeIndex,DynOpt.Ts,DynOpt.OptXstory(:,DynOpt.ActualTimeIndex-1),params);
        DynOpt.OptXstory(:,k) = DynOpt.X; 
        DynOpt.OptXstory_runtime(:,k) = DynOpt.X;

        [temp_Xstory, params] = DynOpt.model_propagate(DynOpt.ActualTimeIndex,DynOpt.Ts,DynOpt.Xstory(:,DynOpt.ActualTimeIndex-1),params);
        DynOpt.Xstory(:,k) = temp_Xstory; 
    else
        DynOpt.OptXstory(:,k) = DynOpt.X_init;
    end

    %%%%%%%%%%%%%%%%%%% MEASUREMENT %%%%%%%%%%%%%%%%%%%%%%%
    % read measure 
    measure_forward = 1;
    [DynOpt.buf_dy,Y_true] = DynOpt.get_measure(DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dy,DynOpt.intY_full_story,params,DynOpt.ActualTimeIndex);
    % copy to Y noise and corrupt only the measure 
    Y_noise = noise_model_v1(Y_true,DynOpt,params);     

    % no filtering
    Y_filter = Y_noise;

    % store total memory
    DynOpt.Ytrue_full_story(:,end+1) = Y_true(:,1);
    DynOpt.Y_full_story(:,end+1) = Y_filter(:,1);
    DynOpt.dY_full_story(:,end+1) = Y_filter(:,2);
    DynOpt.intY_full_story(:,end+1) = Y_filter(:,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k >1
        z = DynOpt.Y_full_story(:,end);
        u = DynOpt.input_true(:,DynOpt.ActualTimeIndex);
        [xnew, Pnew, DynOpt] = Observer_EKF_v2(DynOpt, params,DynOpt.OptXstory(:,DynOpt.ActualTimeIndex-1),z,u);
        DynOpt.OptXstory(:,DynOpt.ActualTimeIndex) = xnew;
        DynOpt.P = Pnew;
    else
        %%% estimated measure %%%
        measure_forward = 1;
        [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params,DynOpt.ActualTimeIndex);
        DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
        DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
        DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%% PERFORMANCE EVALUATION %%%%%%%%%%%%%%%%%
    % fisrt bunch of data - read Y every Nts and check if the signal is
    dJ_cond_v2(DynOpt.theta,DynOpt.beta,DynOpt.gamma);
    
    % clean 
    clc    
    % Display iteration slengthtep
    disp(['Iteration Number: ', num2str(k),'/',num2str(length(DynOpt.time))])
    disp(['Last DJcond: ', num2str(DynOpt.dJ_cond)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end
DynOpt.run_time = toc(run_time);

% output
DynOpt_out = DynOpt;
params_out = params;

end