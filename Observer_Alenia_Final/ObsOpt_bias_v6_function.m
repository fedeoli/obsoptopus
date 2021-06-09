%% ObsOpt implementation
function [DynOpt_out,params_out] = ObsOpt_bias_v6_function(DynOpt,params)

if DynOpt.print
    disp('Processing data with the optimization-based observer...')
end
run_time = tic;
for k=1:length(DynOpt.time)
    
    %%% only iterations
    if k > 1
        fprintf(repmat('\b',1,lineLength))
    end
    lineLength = fprintf('Iteration Number: %d/%d\n',k,length(DynOpt.time));
    %%%%%%%%%%%%%%%%%%

    % update actual index
    DynOpt.ActualTimeIndex = DynOpt.time(k)+1;
    if k == 1 
        if DynOpt.RL_flag == 1
            DynOpt.past_length = size(DynOpt.OptXstory,2)-1;
        else
            DynOpt.past_length = 0;
        end
        if ~isfield(DynOpt,'time_tot')
            DynOpt.time_tot = DynOpt.time;
        end
    end
    DynOpt.Niter = DynOpt.past_length + k;
    DynOpt.Niter_inner = k;

    % reference state - used for noise
    DynOpt.Xtrue = DynOpt.OptXstoryTRUE(:,DynOpt.Niter);

    %forward propagation of the previous estimate
    if(k>1)
        % INTEGRATION OF BOTH POSITION AND ATTITUDE - STACK 
        % Control allocation inside "params" structure          
        [DynOpt.X, params] = DynOpt.model_propagate(DynOpt,DynOpt.Niter,DynOpt.Ts,DynOpt.OptXstory(:,DynOpt.Niter-1),params);
        DynOpt.OptXstory(:,DynOpt.Niter) = DynOpt.X; 
        DynOpt.OptXstory_runtime(:,DynOpt.Niter) = DynOpt.X;

        [temp_Xstory, params] = DynOpt.model_propagate(DynOpt,DynOpt.Niter,DynOpt.Ts,DynOpt.Xstory(:,DynOpt.Niter-1),params);
        DynOpt.Xstory(:,DynOpt.Niter) = temp_Xstory; 
    end

    %%%%%%%%% MEASUREMENT %%%%%%%%%%%
    % read measure 
    measure_forward = 1;

    if DynOpt.simulationModel == 1
        [DynOpt.buf_dY,Y_true] = DynOpt.get_measure(DynOpt,DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dY,DynOpt.intY_full_story,params,DynOpt.Niter);
        % copy to Y noise and corrupt only the measure 
        Y_noise = noise_model_v1(Y_true,DynOpt,params);     
    else
        Y_true = EvaluateCostFunctionOnWindow_Output_general_data(k,measure_forward);
        Y_noise = Y_true;
    end

    % no filtering
    Y_filter = Y_noise;

    % store total memory
    DynOpt.Ytrue_full_story(:,end+1) = Y_true(:,1);
    DynOpt.Y_full_story(:,end+1) = Y_filter(:,1);
    DynOpt.dY_full_story(:,end+1) = Y_filter(:,2);
    DynOpt.intY_full_story(:,end+1) = Y_filter(:,3);

    % fisrt bunch of data - read Y every Nts and check if the signal is
    DynOpt = dJ_cond_v5_function(DynOpt,params);
    distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);
    DynOpt.distance_safe_flag = (distance < DynOpt.safety_interval);
    %%%% select optimisation with hystheresis %%%%%
    hyst_low = (DynOpt.dJ_cond_story(end,max(1,DynOpt.Niter_inner-1)) < DynOpt.dJ_1) && (DynOpt.dJ_cond >= DynOpt.dJ_1);
    hyst_high = (DynOpt.dJ_cond >= DynOpt.dJ_2);
    DynOpt.hyst_flag = ~(hyst_low || hyst_high);
    if  ((distance < DynOpt.Nts) || DynOpt.hyst_flag) && DynOpt.distance_safe_flag
        %%%% ESTIMATED measurements
        % measures
        [DynOpt.buf_dYhat, Yhat] = DynOpt.get_measure(DynOpt,DynOpt.OptXstory(:,DynOpt.Niter),0,measure_forward,DynOpt.buf_dYhat,DynOpt.intYhat_full_story,params,DynOpt.Niter);
        DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
        DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
        DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);

        % clean 
        if DynOpt.print
            clc
        end
    else
        
        if DynOpt.print
            % Display iteration slengthtep
            disp(['n window: ', num2str(DynOpt.w),'  n samples: ', num2str(DynOpt.Nts)])
            disp(['Last cost function: ', num2str(DynOpt.Jstory(end))]);
            disp(['Last DJcond: ', num2str(DynOpt.dJ_cond)]);
            disp(['N. optimisations RUN: ',num2str(DynOpt.opt_counter)]);
            disp(['N. optimisations SELECTED: ',num2str(DynOpt.select_counter)]);
        end
        

        %%%% OUTPUT measurements - buffer of w elements
        % measures
        DynOpt.Y(:,1:end-1) = DynOpt.Y(:,2:end);
        DynOpt.Y(:,end) = Y_filter(:,1);

        % measures derivative
        DynOpt.dY(:,1:end-1) = DynOpt.dY(:,2:end);
        DynOpt.dY(:,end) = Y_filter(:,2);

        % measures integral
        DynOpt.intY(:,1:end-1) = DynOpt.intY(:,2:end);
        DynOpt.intY(:,end) = Y_filter(:,3);

        % backup
        Y_space_backup = DynOpt.Y_space;
        Y_space_full_story_backup = DynOpt.Y_space_full_story;

        % adaptive sampling
        DynOpt.Y_space(1:end-1) = DynOpt.Y_space(2:end);
        DynOpt.Y_space(end) = DynOpt.ActualTimeIndex;
        DynOpt.Y_space_full_story(end+1) = DynOpt.ActualTimeIndex;

        % store measure times
        DynOpt.temp_time = [DynOpt.temp_time k];

        % measures
        [DynOpt.buf_dYhat, Yhat] = DynOpt.get_measure(DynOpt,DynOpt.OptXstory(:,DynOpt.Niter),0,measure_forward,DynOpt.buf_dYhat,DynOpt.intYhat_full_story,params,DynOpt.Niter);
        DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
        DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
        DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %%% forward optimization %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(DynOpt.ForwardOptimization == 1) 

            %%%% FLUSH THE BUFFER IF SAFETY FLAG %%%%
            first_nonzero = find(DynOpt.Y_space,1,'first');
            Y_space_nonzero = DynOpt.Y_space(first_nonzero:end);
            max_dist = max(diff(Y_space_nonzero));
            if isempty(max_dist)
                max_dist = 1;
            end
            n_samples = min(length(DynOpt.Y_space_full_story)-1,DynOpt.w);
            buf_Y_space_full_story = DynOpt.Y_space_full_story(end-n_samples:end);

            % back time index
            buf_dist = diff(buf_Y_space_full_story);
            DynOpt.BackTimeIndex = DynOpt.ActualTimeIndex-sum(buf_dist); 
            DynOpt.BackIterIndex = find(DynOpt.time_tot==DynOpt.BackTimeIndex);

            % set of initial conditions
            if DynOpt.optimise_params == 1
                DynOpt.temp_x0 = DynOpt.OptXstory(DynOpt.integration_pos*6+1:end,DynOpt.BackIterIndex);
            else
                DynOpt.temp_x0 = DynOpt.OptXstory(DynOpt.integration_pos*6+1:end-DynOpt.nparams,DynOpt.BackIterIndex);
            end

            % Optimisation
            % save current state for recovery
            recovery_pos = params.SatellitesCoordinates;
            recovery_att = params.SatellitesAttitude;

            % Optimisation (only if distance_safe_flag == 1)
            opt_time = tic;
            if DynOpt.distance_safe_flag == 1 || DynOpt.always_opt == 1

                % save J before the optimisation to confront it later
                [J_before, ~, DynOpt] = DynOpt.cost_function(DynOpt.temp_x0,params,DynOpt);

                % local params
                params_local = params;

                %%%%% OPTIMISATION - NORMAL MODE %%%%%%
                if DynOpt.fcon_flag == 0
                   [NewXopt, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local,DynOpt),DynOpt.temp_x0,DynOpt.myoptioptions);
                else
                    [NewXopt, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local,DynOpt),DynOpt.temp_x0,DynOpt.Acon, DynOpt.Bcon, DynOpt.Acon_eq, DynOpt.Bcon_eq, DynOpt.lb, DynOpt.ub, ...
                                                    @(x)DynOpt.nonlcon(x,DynOpt,params_local), DynOpt.myoptioptions);
                end

                % use last opt point to update Dynopt buffers
                [~, ~, DynOpt] = DynOpt.cost_function(NewXopt,params,DynOpt);
                
                % opt counter
                DynOpt.opt_counter = DynOpt.opt_counter + 1;
            else
                params.SatellitesCoordinates = recovery_pos;
                params.SatellitesAttitude = recovery_att;
                NewXopt = DynOpt.temp_x0;

                % no care condition 
                J_before = DynOpt.Jstory(end);
                J = J_before;

                % store J chunks
                DynOpt.J_meas_buf = 1;
                DynOpt.J_der_buf = 1;
                DynOpt.J_int_buf = 1;
                DynOpt.J_quat_buf = 1;
                DynOpt.J_spr_buf = 1;
            end

            % adaptive buffer backup restore
            DynOpt.Y_space = Y_space_backup;
            DynOpt.Y_space_full_story = Y_space_full_story_backup;

            % check J_dot condition
            J_diff = (J/J_before);
            distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);

            if ( (J_diff <= DynOpt.Jdot_thresh) || (distance > DynOpt.safety_interval) )  || DynOpt.blue_flag

                % assign optimised state
                if DynOpt.optimise_params == 1
                    DynOpt.X = NewXopt;
                else
                    DynOpt.X = [NewXopt;DynOpt.OptXstory(DynOpt.StateDim+1:end,DynOpt.BackIterIndex)];
                end

                if DynOpt.integration_pos == 1
                    DynOpt.X = [DynOpt.OptXstory(1:6,DynOpt.BackIterIndex); DynOpt.X];
                end

                % store measure times
                DynOpt.opt_chosen_time = [DynOpt.opt_chosen_time k];

                % store J chunks
                DynOpt.J_meas(:,end+1) = DynOpt.J_meas_buf;
                DynOpt.J_der(:,end+1) = DynOpt.J_der_buf;
                DynOpt.J_int(:,end+1) = DynOpt.J_int_buf;
                DynOpt.J_quat(:,end+1) = DynOpt.J_quat_buf;
                DynOpt.J_spr(:,end+1) = DynOpt.J_spr_buf;

                % counters
                DynOpt.jump_flag = 0;
                DynOpt.select_counter = DynOpt.select_counter + 1;

                DynOpt.OptXstory(:,DynOpt.BackIterIndex) = DynOpt.X;

                % params and state update
                if DynOpt.identify == 1
                    [params,DynOpt] = DynOpt.params_update(DynOpt.X,params,DynOpt);
                end
                x_propagate = DynOpt.X;

                %%%%%%%%%%%%%%%%% FIRST MEASURE UPDATE %%%%%%%%
                % manage measurements
                % set the derivative buffer as before the optimisation process (multiple f computation)
                back_time = DynOpt.BackIterIndex;
                if (back_time) >= DynOpt.d1_derivative
                    DynOpt.buf_dYhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
                else
                    init_pos = DynOpt.d1_derivative-back_time;
                    DynOpt.buf_dYhat_temp = [zeros(9,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
                end

                %%%% ESTIMATED measurements
                % measures       
                % NB: the output storage has to be done in
                % back_time+1 as the propagation has been
                % performed 
                [DynOpt.buf_dYhat_temp, Yhat] = DynOpt.get_measure(DynOpt,x_propagate,0,measure_forward,DynOpt.buf_dYhat_temp,DynOpt.intYhat_full_story,params,DynOpt.BackIterIndex);
                DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
                DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
                DynOpt.intYhat_full_story(:,back_time+1) = Yhat(:,3);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

                %%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%
                n_iter_propagate = DynOpt.Niter-DynOpt.BackIterIndex;

                for j =1:n_iter_propagate 
                    % back time
                    back_time = DynOpt.BackIterIndex+j;

                    % integrate
                    [x_propagate, params] = DynOpt.model_propagate(DynOpt,back_time,DynOpt.Ts,x_propagate, params);                      
                    DynOpt.OptXstory(:,back_time) = x_propagate;

                    % manage measurements
                    % set the derivative buffer as before the optimisation process (multiple f computation)
                    if (back_time) >= DynOpt.d1_derivative
                        DynOpt.buf_dYhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
                    else
                        init_pos = DynOpt.d1_derivative-back_time;
                        DynOpt.buf_dYhat_temp = [zeros(9,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
                    end

                    %%%% ESTIMATED measurements
                    % measures       
                    % NB: the output storage has to be done in
                    % back_time+1 as the propagation has been
                    % performed 
                    [DynOpt.buf_dYhat_temp, Yhat] = DynOpt.get_measure(DynOpt,x_propagate,0,measure_forward,DynOpt.buf_dYhat_temp,DynOpt.intYhat_full_story,params,back_time);
                    DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
                    DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
                    DynOpt.intYhat_full_story(:,back_time+1) = Yhat(:,3);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                DynOpt.Jstory(1,end+1) = J;
            else
                % keep the initial guess
                DynOpt.X = DynOpt.temp_x0;
            end                        

            % adaptive sampling
            DynOpt.Y_space(1:end-1) = DynOpt.Y_space(2:end);
            DynOpt.Y_space(end) = DynOpt.ActualTimeIndex;
            DynOpt.Y_space_full_story(end+1) = DynOpt.ActualTimeIndex;

            % stop time counter
            DynOpt.opt_time(end+1) = toc(opt_time);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %%% backward optimization %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else 
            %%%% TO BE DONE %%%%
        end
        
        if DynOpt.print
            clc;
        end
    end
end
DynOpt.run_time = toc(run_time);

% output
DynOpt_out = DynOpt;
params_out = params;

end