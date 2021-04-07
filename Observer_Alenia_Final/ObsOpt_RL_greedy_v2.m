%% greedy algorithm for trajectory optimisation
function RL = ObsOpt_RL_greedy_v2(Nsearch)

    %%%%%%%%%%% INIT SECTION %%%%%%%%%%
    %%% global vars %%%
    global RL struct

    % initialise seed to make results repeatible
    rng(0,'twister');
    
    %%% Define Environment domain %%
    RL.E.domain_target = [pi/4*ones(3,1), pi/4*ones(3,1)];
    RL.E.domain_status = [0.5.*RL.E.domain_target(:,2), 1.5.*RL.E.domain_target(:,2); -deg2rad(10)*ones(3,1), deg2rad(10)*ones(3,1)];
    
    %%% Orbit generation data %%%
    RL.E.domain_ecc = [1e-4; 1.5e-4];
    RL.E.domain_i = [0;pi/2];
    RL.E.domain_om = [0;pi/2];
    RL.E.domain_RAAN = [0;2*pi];
    RL.E.domain_f0 = [0;2*pi];
    RL.E.domain_T = [5.5e3;6.5e3];

    %%%%%%%%%%%%%%%%%%%%%%%%%% ACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Amplitude %%%%
    RL.A.domain_amp = [0, pi/3];
    RL.A.amp_ts = 1e-1;
    RL.A.domain_amp_grid = RL.A.domain_amp(1):RL.A.amp_ts:RL.A.domain_amp(2);
    RL.A.dimAmp = length(RL.A.domain_amp_grid);
    
    %%%% duty cycle %%%%
    RL.A.domain_d = [0, 1];
    RL.A.d_ts = 5e-2;
    RL.A.domain_d_grid = RL.A.domain_d(1):RL.A.d_ts:RL.A.domain_d(2);
    RL.A.dimDuty = length(RL.A.domain_d_grid);
    
    %%%% Period %%%%
    RL.A.domain_T = [10, 100];
    RL.A.T_ts = 10;
    RL.A.domain_T_grid = RL.A.domain_T(1):RL.A.T_ts:RL.A.domain_T(2);
    RL.A.dimT = length(RL.A.domain_T_grid);
    
    %%%% w - buffer %%%%
    RL.A.domain_w = [5, 10];
    RL.A.w_ts = 1;
    RL.A.domain_w_grid = RL.A.domain_w(1):RL.A.w_ts:RL.A.domain_w(2);
    RL.A.dimw = length(RL.A.domain_w_grid);
    
    %%%% Nts %%%%
    RL.A.domain_Nts = [3, 6];
    RL.A.Nts_ts = 1;
    RL.A.domain_Nts_grid = RL.A.domain_Nts(1):RL.A.Nts_ts:RL.A.domain_Nts(2);
    RL.A.dimNts = length(RL.A.domain_Nts_grid);
    
    %%%% Magnetometers %%%%
    RL.A.domain_Magneto = [0, 2];
    RL.A.Magneto_ts = 1;
    RL.A.domain_Magneto_grid = RL.A.domain_Magneto(1):RL.A.Magneto_ts:RL.A.domain_Magneto(2);
    RL.A.dimMagneto = length(RL.A.domain_Magneto_grid);
    
    %%%% Acceptance threshold %%%%
    RL.A.domain_thresh = [0.9, 0.9];
    RL.A.thresh_ts = 0.5;
    RL.A.domain_thresh_grid = RL.A.domain_thresh(1):RL.A.thresh_ts:RL.A.domain_thresh(2);
    RL.A.dimthresh = length(RL.A.domain_thresh_grid);
       
    %%% number of actions field %%%
    RL.A.nActions = 7;
    RL.A.dimActions = RL.A.dimAmp*RL.A.dimDuty*RL.A.dimT*RL.A.dimw*RL.A.dimNts*RL.A.dimMagneto*RL.A.dimthresh;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    init_struct
    struct.RL = 1;

    %%% Greedy initialisation %%%
    % greedy epsilon
    RL.S.eps = 5e-2;

    % Reward estimate init
    RL.S.Q = zeros(RL.A.dimActions,1);

    % Actions counter
    RL.S.N = zeros(RL.A.dimActions,1);
    
    % reward storage
    RL.S.R = zeros(RL.A.dimActions,1);

    % initial guess
    RL.S.bestAction = [RL.A.domain_amp_grid(1); RL.A.domain_d_grid(1); RL.A.domain_T_grid(1); RL.A.domain_w_grid(1); RL.A.domain_Nts_grid(1); RL.A.domain_Magneto_grid(1); RL.A.domain_thresh_grid(1)];
    choice = 'greedy';
    current_reward = 0;

    %%%%%%%%%%% SEARCH PROCEDURE %%%%%%%%%%%%%%
    %%%%%%%%%%%% GREEDY ALGORITHM %%%%%%%%%%%%%
    % search loop
    for i=1:RL.S.Nsearch
        
        % store iteration
        RL.S.i = i;

        %%%%%%%%%%%% random environment definition %%%%%%%%%%%%%%%%%%%%%%%%
        %%% ORBIT %%%
        ecc = (RL.E.domain_ecc(2)-RL.E.domain_ecc(1)).*rand(1) + RL.E.domain_ecc(1);
        inclination = (RL.E.domain_i(2)-RL.E.domain_i(1)).*rand(1) + RL.E.domain_i(1);
        om = (RL.E.domain_om(2)-RL.E.domain_om(1)).*rand(1) + RL.E.domain_om(1);
        RAAN = (RL.E.domain_RAAN(2)-RL.E.domain_RAAN(1)).*rand(1) + RL.E.domain_RAAN(1);
        f0 = (RL.E.domain_f0(2)-RL.E.domain_f0(1)).*rand(1) + RL.E.domain_f0(1);
        T = (RL.E.domain_T(2)-RL.E.domain_T(1)).*rand(1) + RL.E.domain_T(1);
        RL.S.orbit(:,i) = [ecc; inclination; om; RAAN; f0; T];
        
        %%% ATTITUDE %%%
        RL.S.S0(:,i) = (RL.E.domain_status(:,2)-RL.E.domain_status(:,1)).*rand(RL.E.dimState,1) + RL.E.domain_status(:,1);

        %%% select either greedy or explore
        RL.S.p(i) = rand(1);
        if RL.S.p(i) <= RL.S.eps
            % greedy
            RL.S.A(:,i) = RL.S.bestAction;
            choice = 'greedy';
        else
            %%%%%%%% random action choice %%%%%%%%
            %%% amplitude %%%
            pos_amp = randi(RL.A.dimAmp,1);
            temp_amp = RL.A.domain_amp_grid(pos_amp);
            
            %%% duty %%%
            pos_duty = randi(RL.A.dimDuty,1);
            temp_duty = RL.A.domain_d_grid(randi(pos_duty));
            
            %%% period %%%
            pos_T = randi(RL.A.dimT,1);
            temp_T = RL.A.domain_T_grid(pos_T);
            
            %%% w - buffer %%%
            pos_w = randi(RL.A.dimw,1);
            temp_w = RL.A.domain_w_grid(pos_w);
            
            %%% Nts %%%
            pos_Nts = randi(RL.A.dimNts,1);
            temp_Nts = RL.A.domain_Nts_grid(pos_Nts);
            
            %%% Magnetometers %%%
            pos_Magneto = randi(RL.A.dimMagneto,1);
            temp_Magneto = RL.A.domain_Magneto_grid(pos_Magneto);
            
            %%% thresh %%%
            pos_thresh = randi(RL.A.dimthresh,1);
            temp_thresh = RL.A.domain_thresh_grid(pos_thresh);
                        
            RL.S.A(:,i) = [temp_amp; temp_duty; temp_T; temp_w; temp_Nts; temp_Magneto; temp_thresh];
            choice = 'explore';
        end
        
        %%% update action storage %%%
        RL.S.current_pos(i) = pos_amp*RL.A.dimAmp + pos_duty*RL.A.dimDuty + pos_T*RL.A.dimT + pos_w*RL.A.dimw + pos_Nts*RL.A.dimNts + pos_Magneto*RL.A.dimMagneto + pos_thresh*RL.A.dimthresh;
        %%% store current action in memory array %%%
        RL.A.ActionsArray(:,RL.S.current_pos(i)) = RL.S.A(:,i);
        
        %%% display info %%%
        clc
        disp(['Iteration: ', num2str(i), '/', num2str(RL.S.Nsearch)]);
        disp('Available acts: Au           Ad          T             w            Nts          Mag       thresh');
        disp(['Current acts: ', num2str(transpose(RL.S.A(:,i)))]);
        disp(['Current best: ', num2str(transpose(RL.S.bestAction))]);
        disp(['Last action: ', choice]);
        disp(['Last reward: ', num2str(current_reward)]);
        disp('        ');

        % launch satellite simulation
        [DynOpt, params] = ObsOpt_v29_fun_params_alenia(struct);

        % store structures
        RL.S.DynOpt(i) = DynOpt;
        RL.S.params(i) = params;   

        %%% store reward %%
        dJ_cond = RL.S.DynOpt(i).dJ_cond_story(end,:);
        dJ_cond_max = max(RL.S.DynOpt(i).dJ_cond_story(end,:));
        dJ_cond_norm = dJ_cond./dJ_cond_max;
        current_reward = 1/trapz(dJ_cond_norm);
        RL.S.R(RL.S.current_pos(i),end+1) = current_reward;
        
        %%% increase action counter %%%
        RL.S.N(RL.S.current_pos(i)) = RL.S.N(RL.S.current_pos(i))+1;

        %%% update action value %%%
        RL.S.Q(RL.S.current_pos(i)) = RL.S.Q(RL.S.current_pos(i)) + (current_reward-RL.S.Q(RL.S.current_pos(i)))/RL.S.N(RL.S.current_pos(i));

        %%% update best action %%%
        [~,pos] = max(RL.S.Q,[],'all','omitnan','linear');
        RL.S.bestAction = RL.A.ActionsArray(:,pos);
    end
    
    
    %%% BEST ACTION SELECTION %%%
    pos = 0;
    for i=1:RL.S.Nsearch
        if (RL.S.A(:,i) == RL.S.bestAction(:))
            pos = [pos, i];
        end
    end
    RL.S.pos_best = pos(2:end);
    
    
    %%%% PLOT SECTION %%%
    if struct.plot == 1
        plot_RL(RL);
    end
    
    
end


