%% greedy algorithm for trajectory optimisation
function RL = ObsOpt_RL_SARSA_v3(Nsearch,Niter)

    %%%%%%%%%%% INIT SECTION %%%%%%%%%%    
    RL.S.Nsearch = Nsearch;
    RL.S.Niter = Niter;

    % initialise seed to make results repeatible
    rng(0,'twister');
    
    %%%%%%%%%%%%%%%%%%%%%%%% ENVIRONMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Define Environment domain - exploiring start %%%
    %%% define target attitude and initial state range
    RL.E.domain_target = 1*[pi/4*ones(3,1), pi/4*ones(3,1)];
    RL.E.domain_status = [0.5.*RL.E.domain_target(:,2), 1.5.*RL.E.domain_target(:,2); -deg2rad(10)*ones(3,1), deg2rad(10)*ones(3,1)];
    %%% state vector length
    RL.E.dimState = size(RL.E.domain_status,1);
    RL.E.dimTarget = size(RL.E.domain_target,1);
    
    %%% Environment state - nu %%%
    RL.E.domain_S_Ts = 1e-1;
    domain_grid_1d = 0:RL.E.domain_S_Ts:1; 
    RL.E.domain_S = [domain_grid_1d; domain_grid_1d];
    RL.E.domain_S_dim = size(RL.E.domain_S);
    
    %%% Orbit generation data %%%
    RL.E.domain_ecc = [1e-4; 1.5e-4];
    RL.E.domain_i = [0;pi/2];
    RL.E.domain_om = [0;pi/2];
    RL.E.domain_RAAN = [0;2*pi];
    RL.E.domain_f0 = [0;2*pi];
    RL.E.domain_T = [5.5e3;6.5e3];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%% ACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%% Magnetometers %%%%
    RL.A.domain_Magneto = [0, 2];
    RL.A.Magneto_ts = 1;
    RL.A.domain_Magneto_grid = RL.A.domain_Magneto(1):RL.A.Magneto_ts:RL.A.domain_Magneto(2);
    RL.A.dimMagneto = length(RL.A.domain_Magneto_grid);
       
    %%% number of actions field %%%
    RL.A.nActions = 1;
    RL.A.dimActions = RL.A.dimMagneto;  
    RL.A.domain_A = RL.A.domain_Magneto_grid;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %%%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % check row stochastic
    RL.S.pi0 = rand(RL.E.domain_S_dim(1),RL.E.domain_S_dim(2),RL.A.dimActions);
    RL.S.pi0 = RL.S.pi0./sum(RL.S.pi0,3);
    RL.S.Q = 5e-2*RL.S.pi0;
    RL.S.search{1}.pi = RL.S.pi0;
      
    % greedy choice
    RL.S.epsilon = 5e-2;
    
    % reward storage
    RL.S.R = [];
    current_reward = 0;
    
    % algorithm params
    RL.S.alpha = 0.1;
    RL.S.gamma = 0.5;
    
    %%%%%%%%%%% SEARCH PROCEDURE %%%%%%%%%%%%%%
    %%%%%%%%%%%% SARSA ALGORITHM  %%%%%%%%%%%%%
    % search loop
    for s=1:RL.S.Nsearch
        
        %%%%%%%%%%%%%%%%%%%%%%%% INIT EPISODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        init_episode;
        % short period setup
        setup.T_duration = 20;
        setup.t_start = 0;
        setup.Tend = setup.t_start + 20;
        setup.T = 50;
        setup.u_freq = 1/setup.T;
        % init state
        [RL.S.state_pos, RL.S.state] = locate_state(nu,RL.E.domain_S);
        % save setup
        RL.S.setup = setup;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % last search is testing 
        if s == RL.S.Nsearch
            RL.S.epsilon = 1;
            testing = 1;
        else
            testing = 0;
        end
                    
        %%% select action from policy - A
        RL.S.p_eps = rand(1);
        RL.S.p(1) = rand(1);
        if RL.S.p_eps < RL.S.epsilon
            RL.S.action_pos = find(max(getval_user(RL.S.Q,RL.S.state_pos)),1,'first');
            RL.S.A(:,1) = RL.A.domain_A(RL.S.action_pos);
        else
            cum_prob = cumsum(getval_user(RL.S.search{s}.pi,RL.S.state_pos));
            RL.S.action_pos = find(cum_prob>=RL.S.p(1),1,'first');
            RL.S.A(:,1) = RL.A.domain_A(RL.S.action_pos);
        end
        
        % single episode
        for i=1:RL.S.Niter   
            
            % store iteration
            RL.S.i = i;

            %%% display info %%%
            clc
            disp(['Iteration search: ', num2str(s), '/', num2str(RL.S.Nsearch)]);
            disp(['Iteration episode: ', num2str(i), '/', num2str(RL.S.Niter)]);
            disp('Available acts:   Mag');
            disp(['Current acts: ', num2str(transpose(RL.S.A(:,i)))]);
            disp(['Current reward: ', num2str(current_reward)]);
            disp('        ');

            % launch satellite simulation
            setup.RL_data = RL;
            [DynOpt, ~] = ObsOpt_RL_v1_fun(setup);
            DynOpt_save.OptXstory_runtime = DynOpt.OptXstory_runtime;
            DynOpt_save.Opt_quat_runtime = DynOpt.Opt_quat_runtime;
            DynOpt_save.True_quat = DynOpt.True_quat;
            DynOpt_save.OptXstoryTRUE = DynOpt.OptXstoryTRUE;
            DynOpt_save.position_state = DynOpt.position_state;
            DynOpt_save.attitude_state = DynOpt.attitude_state;
            DynOpt_save.time = DynOpt.time;

            % store structures
            RL.S.search{s}.DynOpt{i} = DynOpt_save;  

            %%% evaluate reward %%
            dJ_cond = DynOpt.dJ_cond_story(end,:);
            dJ_cond_max = max(DynOpt.dJ_cond_story(end,:));
            dJ_cond_norm = dJ_cond./dJ_cond_max;
            current_reward = 1/trapz(dJ_cond_norm);
            if isnan(current_reward)
                current_reward = 0;
            end
            RL.S.R(i) = current_reward;

            %%% evaluate next state %%%
            % get last error
            nu = get_state(DynOpt);
            % discretize state
            [RL.S.state_next_pos, RL.S.state_next] = locate_state(nu,RL.E.domain_S);

            %%% select action from policy - A'
            RL.S.p_eps = rand(1);
            RL.S.p(i) = rand(1);
            if RL.S.p_eps < RL.S.epsilon
                pos_1 = find(max(getval_user(RL.S.Q,RL.S.state_pos)),1,'first');
                A_1 = RL.A.domain_A(pos_1);
            else
                cum_prob = cumsum(getval_user(RL.S.search{s}.pi,RL.S.state_pos));
                pos_1 = find(cum_prob>=RL.S.p(i),1,'first');
                A_1 = RL.A.domain_A(pos_1);
            end
            action_next_pos = pos_1;
            action_next = A_1;

            %%% TD0 update value function %%%
            if ~testing
                Q_now = getval_user(RL.S.Q,RL.S.state_pos);
                Q_next = getval_user(RL.S.Q,RL.S.state_next_pos);
                temp =  Q_now(RL.S.action_pos) + RL.S.alpha*(current_reward + RL.S.gamma*Q_next(action_next_pos) - Q_now(RL.S.action_pos));
                pos = [RL.S.state_pos; RL.S.action_pos];
                RL.S.Q = setval_user(RL.S.Q,pos,temp);
            end
            
            % update policy from current Q function
            RL.S.search{s}.pi = RL.S.Q./sum(RL.S.Q,3);
            
            %%% update state %%%
            RL.S.state_pos = RL.S.state_next_pos;
            RL.S.state = RL.S.state_next;

            %%% update action %%%
            RL.S.action_pos = action_next_pos;
            RL.S.A(:,i+1) = action_next;

            %%%%%%%%%%%%%%%%% FOR SATELLITES ONLY - NOT GENERAL %%%%%%%%%%%%%%%
            %%% set correct initial condition for next step %%%
            offset = setup.integration_pos*6;
            RL.S.S0 = DynOpt_save.OptXstory_runtime(offset+1:offset+7,end);
            RL.S.satellites_iner_ECI_true = DynOpt_save.position_state(:,end);
            RL.S.satellites_attitude_true = DynOpt_save.attitude_state(:,end);

            % disable initial error on the estmation, we get the last opt state
            setup.init_error_amp = 0*setup.init_error_amp;
            setup.init_param_error_amp = 0*setup.init_param_error_amp;
            
            setup.t_start = setup.Tend;
            setup.Tend = setup.t_start + 20;
        end
        
        % update policy from current Q function
        RL.S.search{s+1}.pi = RL.S.Q./sum(RL.S.Q,3);
        
        %%% show results
        RL.S.Qval = plotQ(RL.S.Q);
        pause(0.1)
        
    end
    
    %%% show results
    RL.S.Qval = plotQ(RL.S.Q);
    pause(0.1)
end


