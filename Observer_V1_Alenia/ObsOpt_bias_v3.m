%% ObsOpt implementation
function ObsOpt_bias_v3

global DynOpt params

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
        if DynOpt.input_tuning == 0
            % Control allocation inside "params" structure          
            [DynOpt.X, params] = DynOpt.model_propagate(DynOpt.ActualTimeIndex,DynOpt.Ts,DynOpt.OptXstory(:,DynOpt.ActualTimeIndex-1),params);
            DynOpt.OptXstory(:,k) = DynOpt.X; 
        else
            % current time index because it uses the true state
            temp_state = DynOpt.OptXstoryTRUE(1:DynOpt.StateDim+DynOpt.nparams,k);
            DynOpt.X = [temp_state; DynOpt.X(end-2:end)];
            DynOpt.OptXstory(:,k) = DynOpt.X; 
        end
    end

    %%%%%%%%% MEASUREMENT %%%%%%%%%%%
    % read measure 
    measure_forward = 1;

    if DynOpt.simulationModel == 1
        [DynOpt.buf_dy,Y_true] = DynOpt.get_measure(DynOpt.Xtrue,0,measure_forward,DynOpt.buf_dy,DynOpt.intY_full_story,params);
        DynOpt.measure_noise = DynOpt.measure_amp*randn(length(params.observed_state),1);
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
    dJ_cond(DynOpt.theta,DynOpt.beta,DynOpt.gamma);
    distance = DynOpt.ActualTimeIndex-DynOpt.Y_space(end);
    DynOpt.distance_safe_flag = (distance < DynOpt.safety_interval);
    % SENSOR SATURATION/ACCURACY
    DynOpt.sensor_stop_flag = (DynOpt.Y_full_story(DynOpt.ActualTimeIndex) < DynOpt.sensor_stop_thresh);
    % if last sample is too near or the condition isn't met then
    % propagate. However, check that the last optimisation has been
    % performed at most WindowSamples times ago
    if  ((distance < DynOpt.Nts) || (DynOpt.dJ_cond < DynOpt.dJcond_thresh)) && DynOpt.distance_safe_flag
        %%%% ESTIMATED measurements
        % measures
        [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params);
        DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
        DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
        DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);

        % clean 
        clc
    else

        % Display iteration slengthtep
        if DynOpt.input_tuning == 0
            disp(['TARGET: bias: ', num2str(DynOpt.OptXstoryTRUE(end,DynOpt.ActualTimeIndex))])
            disp(['INIT: bias: ', num2str(DynOpt.X_init(8))])
            disp(['CURRENT: bias: ', num2str(params.bias)])
        else
            disp(['INIT: Au: ', num2str(DynOpt.X_init(11)), ' d: ', num2str(DynOpt.X_init(12)), ' f: ', num2str(DynOpt.X_init(13))])
            disp(['CURRENT: Au: ', num2str(DynOpt.X(11)), ' d: ', num2str(DynOpt.X(12)), ' f: ', num2str(DynOpt.X(13))])
        end
        disp(['n window: ', num2str(DynOpt.w),'  n samples: ', num2str(DynOpt.Nts)])
        disp(['Iteration Number: ', num2str(k),'/',num2str(length(DynOpt.time))])
        disp(['Last cost function: ', num2str(DynOpt.Jstory(end))]);
        disp(['N. optimisations RUN: ',num2str(DynOpt.opt_counter)]);
        disp(['N. optimisations SELECTED: ',num2str(DynOpt.select_counter)]);

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

        if (k < max(1,DynOpt.WindowSamples)) 
            [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params);
            DynOpt.Yhat_full_story(:,end+1) = Yhat(:,1);
            DynOpt.dYhat_full_story(:,end+1) = Yhat(:,2);
            DynOpt.intYhat_full_story(:,end+1) = Yhat(:,3);
        else    

            % measures
            [DynOpt.buf_dyhat, Yhat] = DynOpt.get_measure(DynOpt.OptXstory(:,DynOpt.ActualTimeIndex),0,measure_forward,DynOpt.buf_dyhat,DynOpt.intYhat_full_story,params);
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
                if (max_dist >= DynOpt.safety_interval) && (DynOpt.flush_buffer == 1)
                    DynOpt.Y(:,1:end-1) = zeros(DynOpt.dim_out,DynOpt.w-1);
                    DynOpt.dY(:,1:end-1) = zeros(DynOpt.dim_out,DynOpt.w-1);
                    DynOpt.intY(:,1:end-1) = zeros(DynOpt.dim_out,DynOpt.w-1);
                    DynOpt.Y_space(:,1:end-2) = zeros(1,DynOpt.w-2);
                    %%% update also the backup %%%
                    Y_space_backup = zeros(1,DynOpt.w);
                    %%% restore Y_space_full_story
                    n_samples = min(length(DynOpt.Y_space_full_story)-1,DynOpt.w-1);
                    temp_Y_space_full_story = DynOpt.Y_space_full_story;
                    temp_Y_space_full_story(end-n_samples:end-1) = zeros(1,n_samples);
                    buf_Y_space_full_story = temp_Y_space_full_story(end-n_samples:end-1);
                else
                    n_samples = min(length(DynOpt.Y_space_full_story)-1,DynOpt.w);
                    buf_Y_space_full_story = DynOpt.Y_space_full_story(end-n_samples:end);
                end

                % back time index
                buf_dist = diff(buf_Y_space_full_story);
                DynOpt.BackTimeIndex = max(1,k-sum(buf_dist)); 

                % set of initial conditions
                if DynOpt.optimise_params == 1
                    DynOpt.temp_x0 = DynOpt.OptXstory(:,DynOpt.BackTimeIndex);
                else
                    DynOpt.temp_x0 = DynOpt.OptXstory(1:DynOpt.StateDim,DynOpt.BackTimeIndex);
                end

                %%% INPUT TUNING - HYBRID INIT %%%
                if DynOpt.input_tuning == 1
                    input_state = DynOpt.temp_x0(end-2:end-1);
                end

                % Optimisation
                % save current state for recovery
                recovery_pos = params.SatellitesCoordinates;
                recovery_att = params.SatellitesAttitude;

                % Optimisation (only if distance_safe_flag == 1)
                opt_time = tic;
                if DynOpt.distance_safe_flag == 1 || DynOpt.always_opt == 1

                    % save J before the optimisation to confront it later
                    if DynOpt.input_tuning == 0
                        J_before = DynOpt.cost_function(DynOpt.temp_x0,params);
                    else
                        J_before = DynOpt.cost_function(input_state,params);
                    end

                    try
                        % local params
                        params_local = params;

                        if DynOpt.input_tuning == 0
                            %%%%% OPTIMISATION - NORMAL MODE %%%%%%
                            if DynOpt.fcon_flag == 0 
                                [NewXopt, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local),DynOpt.temp_x0,DynOpt.myoptioptions);
                            else
                                if DynOpt.multistart == 0 && DynOpt.globalsearch == 0
                                    [NewXopt, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local),DynOpt.temp_x0,...
                                                                    DynOpt.Acon, DynOpt.Bcon, DynOpt.Acon_eq, DynOpt.Bcon_eq,...
                                                                    DynOpt.lb,DynOpt.ub,DynOpt.nonlcon,DynOpt.myoptioptions);
                                elseif DynOpt.multistart == 1 && DynOpt.globalsearch == 0
                                    problem = createOptimProblem('fmincon','objective',@(x)DynOpt.cost_function(x,params_local),'x0',DynOpt.temp_x0,'lb',DynOpt.lb,'ub',DynOpt.ub,'nonlcon',DynOpt.nonlcon,'options',DynOpt.myoptioptions);
                                    gs = MultiStart;
                                    [NewXopt, J, DynOpt.exitflag] = run(gs,problem,DynOpt.nstart); 
                                elseif DynOpt.globalsearch == 1 && DynOpt.multistart == 0
                                    problem = createOptimProblem('fmincon','objective',@(x)DynOpt.cost_function(x,params_local),'x0',DynOpt.temp_x0,'lb',DynOpt.lb,'ub',DynOpt.ub,'nonlcon',DynOpt.nonlcon,'options',DynOpt.myoptioptions);
                                    gs = GlobalSearch;
                                    [NewXopt, J, DynOpt.exitflag] = run(gs,problem); 
                                else
                                    disp('WRONG OPTIMISATION FLAGS')
                                    return
                                end
                            end
                        else
                            %%% OPTIMISATION - INPUT TUNE MODE %%%%
                            if DynOpt.fcon_flag == 0 
                                [DynOpt.theta_u, J, DynOpt.exitflag] = DynOpt.fmin(@(x)DynOpt.cost_function(x,params_local),input_state,myoptioptions);
                                NewXopt = [DynOpt.temp_x0(1:DynOpt.StateDim+DynOpt.nparams); DynOpt.theta_u; DynOpt.temp_x0(end)];
                            else
                                disp('WRONG OPTIMISATION FLAGS')
                                return
                            end
                        end
                    catch
                        params.SatellitesCoordinates = recovery_pos;
                        params.SatellitesAttitude = recovery_att;
                        NewXopt = DynOpt.temp_x0;
                        J = DynOpt.cost_function(NewXopt,params);
                    end

                    % opt counter
                    DynOpt.opt_counter = DynOpt.opt_counter + 1;
                else
                    params.SatellitesCoordinates = recovery_pos;
                    params.SatellitesAttitude = recovery_att;
                    NewXopt = DynOpt.temp_x0;

                    % no care condition 
                    J_before = 1;
                    J = J_before;
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
                        DynOpt.X = [NewXopt;DynOpt.OptXstory(DynOpt.StateDim+1:end,DynOpt.BackTimeIndex)];
                    end

                    % store measure times
                    DynOpt.opt_chosen_time = [DynOpt.opt_chosen_time k];

                    % counters
                    DynOpt.jump_flag = 0;
                    DynOpt.select_counter = DynOpt.select_counter + 1;

                    if (DynOpt.input_tuning == 0) || (DynOpt.optimise_input == 1)
                        DynOpt.OptXstory(:,DynOpt.BackTimeIndex) = DynOpt.X;

                        % params and state update
                        if DynOpt.identify == 1
                            params = DynOpt.params_update(DynOpt.X,params);
                        end
                        x_propagate = DynOpt.X;

                        %%%%%%%%%%%%%%%%% FIRST MEASURE UPDATE %%%%%%%%
                        % manage measurements
                        % set the derivative buffer as before the optimisation process (multiple f computation)
                        back_time = DynOpt.BackTimeIndex;
                        if (back_time) >= DynOpt.d1_derivative
                            DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
                        else
                            init_pos = DynOpt.d1_derivative-back_time;
                            DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
                        end

                        %%%% ESTIMATED measurements
                        % measures       
                        % NB: the output storage has to be done in
                        % back_time+1 as the propagation has been
                        % performed 
                        [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(x_propagate,0,measure_forward,DynOpt.buf_dyhat_temp,DynOpt.intYhat_full_story,params);
                        DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
                        DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
                        DynOpt.intYhat_full_story(:,back_time+1) = Yhat(:,3);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

                        %%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%
                        n_iter_propagate = DynOpt.ActualTimeIndex-DynOpt.BackTimeIndex;
                        for j =1:n_iter_propagate 
                            % back time
                            back_time = DynOpt.BackTimeIndex+j;

                            % integrate
                            [x_propagate, params] = DynOpt.model_propagate(back_time,DynOpt.Ts,x_propagate, params);                      
                            DynOpt.OptXstory(:,back_time) = x_propagate;

                            % manage measurements
                            % set the derivative buffer as before the optimisation process (multiple f computation)
                            if (back_time) >= DynOpt.d1_derivative
                                DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
                            else
                                init_pos = DynOpt.d1_derivative-back_time;
                                DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
                            end

                            %%%% ESTIMATED measurements
                            % measures       
                            % NB: the output storage has to be done in
                            % back_time+1 as the propagation has been
                            % performed 
                            [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(x_propagate,0,measure_forward,DynOpt.buf_dyhat_temp,DynOpt.intYhat_full_story,params);
                            DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
                            DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
                            DynOpt.intYhat_full_story(:,back_time+1) = Yhat(:,3);
                        end
                    else                               
                        % update with TRUE data
                        DynOpt.OptXstory(end-2:end-1,DynOpt.BackTimeIndex) = DynOpt.theta_u;
                        n_iter_propagate = DynOpt.ActualTimeIndex-DynOpt.BackTimeIndex;
                        for j =1:n_iter_propagate 
                            back_time = DynOpt.BackTimeIndex+j;
                            DynOpt.OptXstory(end-2:end-1,back_time) = DynOpt.theta_u;
                        end

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
        end
        clc;
    end
end
DynOpt.run_time = toc(run_time);

% output
DynOpt_out = DynOpt;
params_out = params;

end