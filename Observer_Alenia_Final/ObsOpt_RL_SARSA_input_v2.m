%% greedy algorithm for trajectory optimisation
function RL = ObsOpt_RL_SARSA_input_v2(Nsearch,test_flag, Q0, Niter)

    %%%%%%%%%%% INIT SECTION %%%%%%%%%%    
    RL.S.Nsearch = Nsearch;

    % initialise seed to make results repeatible
    rng(0,'twister');
    
    %%%%%%%%%%%%%%%%%%%%%%%% ENVIRONMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Define Environment domain - exploiring start %%%
    %%% define target attitude and initial state range
    RL.E.domain_target = 0*[pi/4*ones(3,1), pi/4*ones(3,1)];
    % domain 1 - for nonzero targets
    RL.E.domain_status = [-pi/30*ones(3,1), pi/30*ones(3,1); -deg2rad(10)*ones(3,1), deg2rad(10)*ones(3,1)];
    %%% state vector length
    RL.E.dimState = size(RL.E.domain_status,1);
    RL.E.dimTarget = size(RL.E.domain_target,1);
    
    %%% Environment state - nu %%%
    N_elems = 10;
    RL.E.domain_S_Ts = 1e-1;
    
    % automatic
    domain_grid_1 = linspace(0,1e-1,N_elems); 
    domain_grid_2 = linspace(0,5e-3,N_elems);
    % manual
%     domain_grid_1 = [0, 1e-3, 5e-2, 5e-1, 1]; 
%     domain_grid_2 = [0, 5e-4, 1e-3, 5e-2, 5e-1];
    
    RL.E.domain_S = [domain_grid_1; domain_grid_2];
    RL.E.domain_S_dim = (size(RL.E.domain_S,2)-1)^(size(RL.E.domain_S,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    RL.A.domain_Magneto = [0, 1];
    RL.A.domain_Magneto_grid = RL.A.domain_Magneto;
    RL.A.dimMagneto = length(RL.A.domain_Magneto_grid);
    
    %%%% Input - omega %%%%
    RL.A.step = 5e-5;
    RL.A.Aw_sign = 1;
    RL.A.domain_Aw = [0, 1];
    RL.A.domain_Aw_grid = RL.A.domain_Aw;
    RL.A.dimAw = length(RL.A.domain_Aw_grid);
       
    %%% number of actions field %%%
    RL.A.nActions = 2;
    RL.A.dimActions = RL.A.dimAw*RL.A.dimMagneto;  
    RL.A.domain_A = permn(0:RL.A.nActions-1,RL.A.nActions)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %%%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % check row stochastic
    if test_flag == 0
        RL.S.Q = 5e-2*rand(RL.E.domain_S_dim,RL.A.dimActions);
        RL.S.Q(1,:) = [0 0 0 1];
    else
        RL.S.Q = Q0;
    end    
    RL.S.search{1}.Q = RL.S.Q;
    
    % N count
    RL.S.N_state = zeros(size(RL.E.domain_S,2)-1);
      
    % greedy choice 
    if test_flag
        RL.S.epsilon = 1;
    else
        RL.S.epsilon = 0.95;
    end
        
    % reward storage
    RL.S.R = [];
    current_reward = 0;
    
    % algorithm params
    RL.S.alpha = 0.1;
    RL.S.gamma = 0.5;
    
    % first RL assignmement
    setup.RL_data = RL;
    
    % Niter
    if test_flag
       RL.S.Niter = Niter; 
    else
       RL.S.Niter = 0; 
    end
    
    %%%%%%%%%%% SEARCH PROCEDURE %%%%%%%%%%%%%%
    %%%%%%%%%%%% SARSA ALGORITHM  %%%%%%%%%%%%%
    % search loop
    for s=1:RL.S.Nsearch
        
        % store iteration
        RL.S.s = s;
        
        %%%%%%%%%%%%%%%%%%%%%%%% INIT EPISODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        init_flag = 1;
        init_episode;
        init_flag = 0;
        % init state
        [RL.S.state_pos_ind, RL.S.state_pos, RL.S.state] = locate_state(nu,RL.E.domain_S);
        % save setup
        RL.S.setup = setup;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
        %%% select action from policy - A
        RL.S.p(1) = rand(1);
        if RL.S.p(1) < RL.S.epsilon
            RL.S.action_pos = find(getval_user(RL.S.Q,RL.S.state_pos_ind) == max(getval_user(RL.S.Q,RL.S.state_pos_ind)),1,'first');
            RL.S.A(:,1) = RL.A.domain_A(:,RL.S.action_pos);
        else
            RL.S.action_pos = randi(RL.A.dimActions);
            RL.S.A(:,1) = RL.A.domain_A(:,RL.S.action_pos);
        end
        
        % single episode
        RL.S.i = 1;
        terminal = 0;
        [~, sub_pos, ~] = locate_state(nu,RL.E.domain_S);
        
        % SARSA loop
        while (terminal == 0) && (test_flag == 0) || ((RL.S.i < RL.S.Niter) && (test_flag == 1))
            
            % store iteration
            i = RL.S.i;

            %%% display info %%%
            clc
            if test_flag
                disp(['Testing - Iteration: ', num2str(RL.S.i),'/',num2str(RL.S.Niter)]);
            else
                disp(['Learning - Episode: ', num2str(s), '/', num2str(RL.S.Nsearch)]);
            end
            disp(['Current state: ', num2str(sub_pos')]);
            disp(['Current acts: ', num2str(transpose(RL.S.A(:,i)))]);
            disp(['Current reward: ', num2str(current_reward)]);
            disp('        ');

            % launch satellite simulation
            setup.RL_data = RL;
            try 
                [DynOpt, ~] = ObsOpt_RL_v1_fun(setup);
                fail_flag = 0;
            catch ME
                err_str = getReport(ME);   
                fail_flag = 1;
            end
            if sign(DynOpt.Aw(1))
                RL.S.Aw_sign = sign(DynOpt.Aw(1));
            end
            
            if fail_flag == 0
                DynOpt_save.OptXstory_runtime = DynOpt.OptXstory_runtime;
                DynOpt_save.Opt_quat_runtime = DynOpt.Opt_quat_runtime;
                DynOpt_save.True_quat = DynOpt.True_quat;
                DynOpt_save.OptXstoryTRUE = DynOpt.OptXstoryTRUE;
                DynOpt_save.position_state = DynOpt.position_state;
                DynOpt_save.attitude_state = DynOpt.attitude_state;
                DynOpt_save.time = DynOpt.time;
                DynOpt_save.OptErrorStory_Euler = DynOpt.OptErrorStory_Euler;
                DynOpt_save.desatt_true = DynOpt.desatt_true;
            else
                DynOpt_save = 'failed run';
                DynOpt.dJ_cond = Inf;
                DynOpt.dJ_cond_story = Inf*ones(4,2);
                DynOpt.track_err = nan;                
            end

            % store structures
            RL.S.search{s}.DynOpt{i} = DynOpt_save;  

            %%% evaluate reward %%
            nu_R = get_state_v2(DynOpt);
            [ind_pos, sub_pos, ~] = locate_state(nu_R,RL.E.domain_S);
            if (sub_pos(1) == 1) && (sub_pos(2) == 1)
                current_reward = 0;
                terminal = 1;
            else
                current_reward = -1;
                terminal = 0;
            end
            RL.S.R(i) = current_reward;
            RL.S.sub_pos(:,i) = sub_pos;
            RL.S.N_state(ind_pos) = RL.S.N_state(ind_pos) + 1;

            %%% evaluate next state %%%
            % get last error
            nu = get_state_v2(DynOpt);
            % discretize state
            [RL.S.state_next_pos_ind, RL.S.state_next_pos, RL.S.state_next] = locate_state(nu,RL.E.domain_S);

            %%% select action from policy - A'
            RL.S.p(i) = rand(1);
            if RL.S.p(i) < RL.S.epsilon
                pos_1 = find(getval_user(RL.S.Q,RL.S.state_next_pos_ind) == max(getval_user(RL.S.Q,RL.S.state_next_pos_ind)),1,'first');
                A_1 = RL.A.domain_A(:,pos_1);
            else
                pos_1 = randi(RL.A.dimActions);
                A_1 = RL.A.domain_A(:,pos_1);
            end
            action_next_pos = pos_1;
            action_next = A_1;

            %%% TD0 update value function %%%
            if ~test_flag
                Q_now = getval_user(RL.S.Q,RL.S.state_pos_ind);
                Q_next = getval_user(RL.S.Q,RL.S.state_next_pos_ind);
                temp =  Q_now(RL.S.action_pos) + RL.S.alpha*(current_reward + RL.S.gamma*Q_next(action_next_pos) - Q_now(RL.S.action_pos));
                pos = [RL.S.state_pos_ind; RL.S.action_pos];
                RL.S.Q = setval_user(RL.S.Q,pos,temp);
            end
            
            %%% update state %%%
            RL.S.state_pos_ind = RL.S.state_next_pos_ind;
            RL.S.state = RL.S.state_next;

            %%% update action %%%
            RL.S.action_pos = action_next_pos;
            RL.S.A(:,i+1) = action_next;
            
            % update iter %
            RL.S.i = RL.S.i + 1;
            
            if fail_flag == 0
                %%% FOR SATELLITES ONLY - NOT GENERAL %%%
                [setup,RL] = save_past_sat(setup,DynOpt,DynOpt_save,RL);  
            else
                terminal = 1;
            end
        end
        
        % update policy from current Q function
        RL.S.search{s+1}.Q = RL.S.Q;
                        
        %%% show results
        RL.S.Qval = plotQ(RL.S.Q,RL.E.domain_S,RL.A.domain_A,RL.S.N_state);
        
    end
    
    %%% show results
    RL.S.Qval = plotQ(RL.S.Q,RL.E.domain_S,RL.A.domain_A,RL.S.N_state);
end


