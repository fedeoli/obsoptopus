%% greedy algorithm for trajectory optimisation
function RL = ObsOpt_RL_SARSA_FA(Nsearch,test_flag, w0, Niter)

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
    RL.A.Aw_sign = ones(3,1);
    RL.A.domain_Aw = [0, 1];
    RL.A.domain_Aw_grid = RL.A.domain_Aw;
    RL.A.dimAw = length(RL.A.domain_Aw_grid);
       
    %%% number of actions field %%%
    RL.A.nActions = 4;
    RL.A.dimActions = RL.A.dimAw*RL.A.dimMagneto;  
    temp = permn(0:RL.A.nActions-1,RL.A.nActions)';
    [ ~ , pos_false] = find(temp(1,:) > 1);
    temp(:,pos_false) = [];
    for i=1:3
        [ ~ , pos_false] = find(temp(i+1,:) > 1);
        temp(:,pos_false) = [];
    end
    RL.A.domain_A = temp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %%%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
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
    
    % function approximation init
    if test_flag
        RL.S.w = w0;
    else
        RL.S.w = rand(7,1);
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
        nu = get_state_v2(DynOpt);
        terminal = isterminal(nu);
        init_flag = 0;
        % save setup
        RL.S.setup = setup;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
        %%% select action from policy - A
        RL.S.p(1) = rand(1);
        if RL.S.p(1) < RL.S.epsilon
            q_val = zeros(1,RL.A.dimActions);
            for a=1:RL.A.dimActions
                q_val(a) = q_fun(nu,a,RL.S.w);
            end
            RL.S.action_pos = find(q_val == max(q_val),1,'first');
            RL.S.A(:,1) = RL.A.domain_A(:,RL.S.action_pos);
        else
            RL.S.action_pos = randi(RL.A.dimActions);
            RL.S.A(:,1) = RL.A.domain_A(:,RL.S.action_pos);
        end
        
        % single episode
        RL.S.i = 1;
        
        % SARSA loop
        while (terminal == 0) && (test_flag == 0) || ((RL.S.i < RL.S.Niter) && (test_flag == 1))
            
            % store iteration
            i = RL.S.i;
            
            % store current state
            nu = get_state_v2(DynOpt);
            RL.S.N_state(:,i) = nu;

            %%% display info %%%
            clc
            if test_flag
                disp(['Testing - Iteration: ', num2str(RL.S.i),'/',num2str(RL.S.Niter)]);
            else
                disp(['Learning - Episode: ', num2str(s), '/', num2str(RL.S.Nsearch)]);
            end
            disp(['Current state: ', string(nu')]);
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
            
            % get input sign
            RL.A.Aw_sign = sign(DynOpt.Aw);
            
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

            %%% evaluate reward and next state %%
            nu_R = get_state_v2(DynOpt);
            terminal = isterminal(nu_R);
            if terminal
                current_reward = 0;
            else
                current_reward = -1;
            end
            RL.S.R(i) = current_reward;

            %%% select action from policy - A'
            RL.S.p(i) = rand(1);
            if RL.S.p(i) < RL.S.epsilon
                q_val = zeros(1,RL.A.dimActions);
                for a=1:RL.A.dimActions
                    q_val(a) = q_fun(nu_R,a,RL.S.w);
                end
                pos_1 = find(q_val == max(q_val),1,'first');
            else
                pos_1 = randi(RL.A.dimActions);
            end
            A_1 = RL.A.domain_A(:,pos_1);
            action_next_pos = pos_1;    
            action_next = A_1;

            %%% TD0 update value function %%%
            if ~test_flag
                try
                RL.S.w = RL.S.w + RL.S.alpha*(current_reward + RL.S.gamma*q_fun(nu_R,action_next_pos,RL.S.w) - q_fun(nu,RL.S.action_pos,RL.S.w))*[nu;RL.S.action_pos];
                catch
                   a=1; 
                end
            end

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
        RL.S.search{s+1}.w = RL.S.w;
                               
    end
end


