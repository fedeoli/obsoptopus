%% init episode for RL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DYNOPT INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initialise DynOpt %%%%     
init_struct
setup.load_mem = 0;
%%%% setup as it was the first algorithm iteration (MDP state init) %%%
setup.RL = 1;
% state duration and main setup
setup.w = 5;
setup.Nts = 3;
setup.d = 0;
setup.u_amp = 0;
% params estimation
setup.optimise_params = 1;
setup.nparams = 1;
setup.inertia = 1;
% iteration
RL.S.i = 1;
% action
RL.S.A(:,1) = RL.A.domain_A(:,1);
% initial attitude - true
RL.S.satellites_attitude_true = (RL.E.domain_status(:,2)-RL.E.domain_status(:,1)).*rand(RL.E.dimState,1) + RL.E.domain_status(:,1);
% initial attitude - est
setup.noise_enable = 1;
% error init definition
setup.ParamsAllPos = 1;
setup.init_param_error_amp = 1*setup.noise_enable*ones(1,setup.nparams)*50;
% init state
if setup.noise_enable
    RL.S.S0 = (RL.E.domain_status(:,2)-RL.E.domain_status(:,1)).*rand(RL.E.dimState,1) + RL.E.domain_status(:,1);
else
    RL.S.S0 = RL.S.satellites_attitude_true;
end 
% target attitude
RL.S.T0 = (RL.E.domain_target(:,2)-RL.E.domain_target(:,1)).*rand(RL.E.dimTarget,1) + RL.E.domain_target(:,1);
setup.RL_data = RL;
%%%% get first nu and iner_ECI %%%%
setup.Tend = 2;
[DynOpt, ~] = ObsOpt_RL_v1_fun(setup);
nu = get_state(DynOpt);

%%%%%%%%%%%% random environment definition %%%%%%%%%%%%%%%%%%%%%%%%
%%% ORBIT %%%
ecc = (RL.E.domain_ecc(2)-RL.E.domain_ecc(1)).*rand(1) + RL.E.domain_ecc(1);
inclination = (RL.E.domain_i(2)-RL.E.domain_i(1)).*rand(1) + RL.E.domain_i(1);
om = (RL.E.domain_om(2)-RL.E.domain_om(1)).*rand(1) + RL.E.domain_om(1);
RAAN = (RL.E.domain_RAAN(2)-RL.E.domain_RAAN(1)).*rand(1) + RL.E.domain_RAAN(1);
f0 = (RL.E.domain_f0(2)-RL.E.domain_f0(1)).*rand(1) + RL.E.domain_f0(1);
T = (RL.E.domain_T(2)-RL.E.domain_T(1)).*rand(1) + RL.E.domain_T(1);
RL.S.orbit(:,RL.S.i) = [ecc; inclination; om; RAAN; f0; T];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%