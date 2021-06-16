%% greedy algorithm for trajectory optimisation
function RL = ObsOpt_RL_TD_FA(Nsearch,lambda,test_flag, w0, Niter)

    %%%%%%%%%%% INIT SECTION %%%%%%%%%%    
    RL.S.Nsearch = Nsearch;

    % initialise seed to make results repeatible
    rng(0,'twister');
    
    %%%%%%%%%%%%%%%%%%%%%%%% ENVIRONMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Define Environment domain - exploiring start %%%
    %%% define target attitude and initial state range
    RL.E.domain_target = 0*[pi/4*ones(3,1), pi/4*ones(3,1)];
    % domain 1 - for nonzero targets
    RL.E.domain_status = 0*[-pi/4*ones(3,1), pi/4*ones(3,1); -deg2rad(10)*ones(3,1), deg2rad(10)*ones(3,1)];
    %%% state vector length
    RL.E.dimState = size(RL.E.domain_status,1);
    RL.E.dimTarget = size(RL.E.domain_target,1);
    
    %%% Orbit generation data %%%
    RL.E.domain_ecc = [1.2720e-04; 1.2720e-04];
    RL.E.domain_i = 1*[98.1829*pi/180; 98.1829*pi/180];
    RL.E.domain_om = 1*[85.7508*pi/180;85.7508*pi/180];
    RL.E.domain_RAAN = 1*[208.3314*pi/180;208.3314*pi/180];
    RL.E.domain_f0 = 1*[0*pi/180;0*pi/180];
    RL.E.domain_T = [5.5e3;6.5e3];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%% ACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%% Input - omega %%%%
    RL.A.step = 2e-5;
    RL.A.Aw_sign = ones(3,1);
    RL.A.domain_Aw = [0, 1];
    RL.A.domain_Aw_grid = RL.A.domain_Aw;
    RL.A.dimAw = length(RL.A.domain_Aw_grid);
    
    %%%% Input - omega %%%%
    RL.A.step_mag = 1e-2;
    RL.A.Amag = 0.5;%RL.A.step_mag*randi(floor(1/RL.A.step_mag));
    RL.A.domain_Amag = [-1, 1];
    RL.A.domain_Amag_grid = RL.A.domain_Amag;
    RL.A.dimAmag = length(RL.A.domain_Amag_grid);
       
    %%% number of actions field %%%
    RL.A.nActions = 1;

    RL.A.domain_A = permn(1:2,RL.A.nActions)';
    RL.A.dimActions = size(RL.A.domain_A,2);  
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
    RL.S.alpha = 1e-2;
    RL.S.gamma = 0.5;
    RL.S.lambda = lambda;
    
    % first RL assignmement
    setup.RL_data = RL;
    
    % Niter
    if test_flag
       RL.S.Niter = Niter; 
    else
       RL.S.Niter = 0; 
    end
    
    % function approximation init
    RL.S.M = 5;
    RL.S.N = 4;
    RL.S.statedim = 2;
    RL.S.terminal_cond = [1e-3*ones(1,1); 1e-3*ones(1,1)];
    RL.S.terminal_streak_cond = 5;
    RL.S.terminal_streak = 0;
    RL.S.tile = init_tile(RL.S.statedim,RL.A.dimActions,RL.S.M,RL.S.N);
    RL.S.d = RL.A.dimActions*RL.S.N*(RL.S.M+1)^RL.S.statedim;
    if test_flag
        RL.S.w = w0;
    else
        RL.S.w = rand(RL.S.d,1);
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
        RL = isterminal(nu,RL.S.terminal_cond,RL.S.terminal_streak_cond,RL);
        init_flag = 0;
        % save setup
        RL.S.setup = setup;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % single episode
        RL.S.i = 1;
        
        % elegibility traces
        RL.S.z = zeros(RL.S.d,1);
        Vold = 0;
        
        % init action
        RL.S.action_pos(1) = randi(RL.A.dimActions);
        RL.S.A(:,1) = RL.A.domain_A(:,RL.S.action_pos(1));
        
        % SARSA loop
        while (RL.S.isTerminal == 0) && (test_flag == 0) || ((RL.S.i <= RL.S.Niter) && (test_flag == 1))
           
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
%             action = [DynOpt.Aw(1), DynOpt.y_weight(4)];
            action = DynOpt.Aw(1);
            disp(['Current acts: ', num2str(action)]);
            disp(['Current reward: ', num2str(current_reward)]);
            disp(['Terminal streak: ', num2str(RL.S.terminal_streak),'/',num2str(RL.S.terminal_streak_cond)])
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
            % get Mag amplitude
            RL.A.Amag = DynOpt.y_weight(4);
            
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
            RL = isterminal(nu_R,RL.S.terminal_cond,RL.S.terminal_streak_cond,RL);
            if RL.S.isTerminal == 1
                if fail_flag == 0
                    current_reward = 0;
                else
                    current_reward = -100;
                end
            else
                current_reward = -1;
            end
            RL.S.R(i) = current_reward;
            
            %%% select action from policy - A'
            RL.S.p(1) = rand(1);
            if RL.S.p(1) < RL.S.epsilon
                q_val = zeros(1,RL.A.dimActions);
                for a=1:RL.A.dimActions
                    [x_greedy,~] = getFeatures(RL.S.tile,RL.S.M,RL.S.N,RL.A.dimActions,a,nu);
                    q_val(a) = transpose(RL.S.w)*x_greedy;
                end
                action_pos = find(q_val == max(q_val),1,'first');
                action_next = RL.A.domain_A(:,action_pos);
            else
                action_pos = randi(RL.A.dimActions);
                action_next = RL.A.domain_A(:,action_pos);
            end
                        
            %%% TD store data %%%
            if ~test_flag
                [x,x_ind] = getFeatures(RL.S.tile,RL.S.M,RL.S.N,RL.A.dimActions,RL.S.action_pos(i),nu);
                [x_R,x_R_ind] = getFeatures(RL.S.tile,RL.S.M,RL.S.N,RL.A.dimActions,action_pos,nu_R);
                V = transpose(RL.S.w)*x;
                V_R = transpose(RL.S.w)*x_R;
                delta = current_reward + RL.S.gamma*V_R - V;
                RL.S.z = RL.S.gamma*RL.S.lambda*RL.S.z + (1-RL.S.alpha*RL.S.gamma*RL.S.lambda*transpose(RL.S.z)*x)*x;
                RL.S.w = RL.S.w + RL.S.alpha*(delta + V - Vold)*RL.S.z - RL.S.alpha*(V - Vold)*x;
                Vold = V_R;
            end
            
            % update iter %
            RL.S.i = RL.S.i + 1;
            
            % update action
            RL.S.action_pos(i+1) = action_pos;
            RL.S.A(:,i+1) = action_next;
            
            if fail_flag == 0
                %%% FOR SATELLITES ONLY - NOT GENERAL %%%
                [setup,RL] = save_past_sat(setup,DynOpt,DynOpt_save,RL);  
            else
                RL.S.isTerminal = 1;
            end
        end
        
        % update policy from current Q function
            RL.S.search{s+1}.w = RL.S.w;          
    end    
end


