%% greedy algorithm for trajectory optimisation
function RL = ObsOpt_RL_greedy_v2(Nsearch)

    %%%%%%%%%%% INIT SECTION %%%%%%%%%%
    %%% global vars %%%
    global RL struct

    % initialise seed to make results repeatible
    rng(0,'twister');
    
    %%% Define Environment domain %%
    RL.E.domain_T = [200, 200];
    RL.E.domain_freq = 1./RL.E.domain_T;
    RL.E.domain_target = [-pi*ones(3,1), pi*ones(3,1)];
    RL.E.domain_status = [0.5.*RL.E.domain_target(:,2), 1.5.*RL.E.domain_target(:,2); -deg2rad(10)*ones(3,1), deg2rad(10)*ones(3,1)];

    %%% Define actions domain %%
    % Amplitude 
    RL.A.domain_amp = [0, pi/2];
    RL.A.amp_ts = 1e-1;
    RL.A.domain_amp_grid = RL.A.domain_amp(1):RL.A.amp_ts:RL.A.domain_amp(2);
    % duty cycle
    RL.A.domain_d = [0, 1];
    RL.A.d_ts = 5e-2;
    RL.A.domain_d_grid = RL.A.domain_d(1):RL.A.d_ts:RL.A.domain_d(2);

    % actions length
    RL.A.dimAmp = length(RL.A.domain_amp_grid);
    RL.A.dimDuty = length(RL.A.domain_d_grid);

    % search length setup
    RL.S.Nsearch = Nsearch;

    % status length
    RL.E.dimState = size(RL.E.domain_status,1);
    RL.E.dimTarget = size(RL.E.domain_target,1);

    % storage
    RL.S.S0 = zeros(RL.E.dimState, RL.S.Nsearch);
    RL.S.F0 = zeros(1,RL.S.Nsearch);
    RL.S.T0 = zeros(RL.E.dimTarget,RL.S.Nsearch);


    %%%% initialise DynOpt %%%%     
    init_struct_RL_v2

    %%% Greedy initialisation %%%
    % greedy epsilon
    RL.S.eps = 5e-2;

    % Reward estimate init
    RL.S.Q = zeros(RL.A.dimAmp, RL.A.dimDuty);

    % Actions counter
    RL.S.N = zeros(RL.A.dimAmp, RL.A.dimDuty);
    
    % reward storage
    RL.S.R = zeros(RL.A.dimAmp, RL.A.dimDuty,1);

    % initial guess
    RL.S.bestAction = [RL.A.domain_amp_grid(1); RL.A.domain_d_grid(1)];
    choice = 'greedy';
    pos_amp = 1;
    pos_duty = 1;
    current_reward = 0;

    %%%%%%%%%%% SEARCH PROCEDURE %%%%%%%%%%%%%%
    %%%%%%%%%%%% GREEDY ALGORITHM %%%%%%%%%%%%%
    % search loop
    for i=1:RL.S.Nsearch
        
        %%% display info %%%
        clc
        disp(['Iteration: ', num2str(i), '/', num2str(RL.S.Nsearch)]);
        disp(['Current best: Amp = ', num2str(RL.S.bestAction(1)), '  Duty = ', num2str(RL.S.bestAction(2))]);
        disp(['Last action: ', choice]);
        disp(['Last reward: ', num2str(current_reward)]);
        disp('        ');
        
        % store iteration
        RL.S.i = i;

        % random environment definition
        RL.S.S0(:,i) = (RL.E.domain_status(:,2)-RL.E.domain_status(:,1)).*rand(RL.E.dimState,1) + RL.E.domain_status(:,1);
        RL.S.F0(i) = (RL.E.domain_freq(2)-RL.E.domain_freq(1)).*rand(1,1) + RL.E.domain_freq(1);
        RL.S.T0(:,i) = (RL.E.domain_target(:,2)-RL.E.domain_target(:,1)).*rand(RL.E.dimTarget,1) + RL.E.domain_target(:,1);

        %%% select either greedy or explore
        RL.S.p(i) = rand(1);
        if RL.S.p(i) <= RL.S.eps
            % greedy
            RL.S.A(:,i) = RL.S.bestAction;
            choice = 'greedy';
        else
            % random action choice
            pos_amp = randi(RL.A.dimAmp,1);
            pos_duty = randi(RL.A.dimDuty,1);
            temp_amp = RL.A.domain_amp_grid(pos_amp);
            temp_duty = RL.A.domain_d_grid(randi(pos_duty));
            RL.S.A(:,i) = [temp_amp; temp_duty];
            choice = 'explore';
        end

        % launch satellite simulation
        [DynOpt, params] = ObsOpt_RL_v3(struct,RL);

        % store structures
        RL.S.DynOpt(i) = DynOpt;
        RL.S.params(i) = params;   

        %%% store reward %%
        current_reward = 1/trapz(RL.S.DynOpt(i).dJ_cond_story(end,:));
        RL.S.R(pos_amp,pos_duty,end+1) = current_reward;
        
        %%% increase action counter %%%
        RL.S.N(pos_amp,pos_duty) = RL.S.N(pos_amp,pos_duty)+1;

        %%% update action value %%%
        RL.S.Q(pos_amp,pos_duty) = RL.S.Q(pos_amp,pos_duty) + (current_reward-RL.S.Q(pos_amp,pos_duty))/RL.S.N(pos_amp,pos_duty);

        %%% update best action %%%
        [~,pos] = max(RL.S.Q,[],'all','omitnan','linear');
        [row, col] = ind2sub(size(RL.S.Q),pos);
        RL.S.bestAction = [RL.A.domain_amp_grid(row); RL.A.domain_d_grid(col)];
    end
    
    %%% POST LEARNING DATA ANALYSYS %%%
    for i=1:RL.A.dimAmp
        for j=1:RL.A.dimDuty
            RL.S.sigma(i,j) = std(nonzeros(RL.S.R(i,j,:)),'omitnan');
        end
    end
    
    %%% BEST ACTION SELECTION %%%
    pos = 0;
    for i=1:RL.S.Nsearch
        if (RL.S.DynOpt(i).u_amp == RL.S.bestAction(1)) && (RL.S.DynOpt(i).d == RL.S.bestAction(2))
            pos = [pos, i];
        end
    end
    RL.S.pos_best = pos(2:end);
    
    %%% COMPARISON ACTION SELECTION %%%
    pos = 0;
    for i=1:RL.S.Nsearch
        if (RL.S.DynOpt(i).u_amp == RL.A.domain_amp_grid(3)) && (RL.S.DynOpt(i).d == RL.A.domain_d_grid(5))
            pos = [pos, i];
        end
    end
    RL.S.pos_compare = pos(2:end);
    
    %%%% PLOT SECTION %%%
    if struct.plot == 1
        plot_RL(RL);
    end
    
    
end


