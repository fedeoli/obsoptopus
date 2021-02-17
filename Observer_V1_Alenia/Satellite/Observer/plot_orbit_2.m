%% plot orbits - from official data
%% new agents position error
if 1
    
work1 = 'simulations/Init_near/4_500s_OldNoise/GPS_EKF_new';
work2 = 'simulations/Init_near/4_500s_OldNoise/GPS_UKF_new';

%%%% KALMAN FILTERS SECTION
load(work1)
traj_error_1 = traj_error_KF;
traj_error_1_sum = traj_error_KF_sum;
mean_traj_error_1_sum = mean_traj_error_KF_sum;
sigma_traj_error_1_sum = sigma_traj_error_KF_sum;

load(work2)
traj_error_2 = traj_error_KF;
traj_error_2_sum = traj_error_KF_sum;
mean_traj_error_2_sum = mean_traj_error_KF_sum;
sigma_traj_error_2_sum = sigma_traj_error_KF_sum;

%%%% GPS SECTION
% load(work1)
% traj_error_1 = traj_error_GPS;
% traj_error_1_sum = traj_error_GPS_sum;
% mean_traj_error_1_sum = mean_traj_error_GPS_sum;
% sigma_traj_error_1_sum = sigma_traj_error_GPS_sum;
% 
% load(work2) 
% traj_error_2 = traj_error_GPS;
% traj_error_2_sum = traj_error_GPS_sum;
% mean_traj_error_2_sum = mean_traj_error_GPS_sum;
% sigma_traj_error_2_sum = sigma_traj_error_GPS_sum;

%%%% 3rd PLOT SECTION
% load simulations\4_500s\GPS
% traj_error_3 = traj_error_GPS;
% traj_error_3_sum = traj_error_GPS_sum;
% mean_traj_error_3_sum = mean_traj_error_GPS_sum;
% sigma_traj_error_3_sum = sigma_traj_error_GPS_sum;

nagent = ObserverTest.Nagents;
% start_step_local = ObserverTest.window_interval(1);
start_step_local = 60;

% Observerplot plots
if 1
    figure('Name', 'Agents position estimation errors (norm)')
    hold on
    ylabel('Position: mean+-sigma [m]')
    xlabel('Agents')
    grid on 
    
    load(work1)
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.Mean(ObserverTest.Ntest,:)*1000,ObserverTest.Sigma(ObserverTest.Ntest,:)*1000,'LineWidth',2);

    load(work2)
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.Mean(ObserverTest.Ntest,:)*1000,ObserverTest.Sigma(ObserverTest.Ntest,:)*1000,'LineWidth',2);
%     load simulations/Init_near/GPSonly/GPS_HN_6000
%     errorbar(1:1:ObserverTest.Nagents, ObserverTest.Mean(ObserverTest.Ntest,:)*1000,ObserverTest.Sigma(ObserverTest.Ntest,:)*1000,'LineWidth',2);
    
%     legend('EKF','GPS')
%     legend('EKF','GPS+EKF')
%     legend('UKF','GPS+UKF')
%     legend('EKF','UKF','GPS')
%     legend('GPS+EKF','GPS+UKF')
    legend('GPS+EKF+UWB','GPS+UKF+UWB')
%     legend('GPS','GPS opt')
%     legend('GPS_UKF','GPS')
%     legend('UKF','GPS+UKF')
%     legend('GPS')
%     legend('Chi','GPS')
%     legend('GPS+UKF','GPS+UKF 8 ag')
end

if 1
    figure
    ylabel('Position: norm of agents error [Km]')
    xlabel('Integration steps')
    hold on
    grid on
    
    plot(start_step_local:ObserverTest.Npassi,traj_error_1_sum(start_step_local:end) ,'b-','LineWidth',2);
    plot(start_step_local:ObserverTest.Npassi,traj_error_2_sum(start_step_local:end) ,'r-','LineWidth',2);
%     plot(start_step:ObserverTest.Npassi,traj_error_3_sum(start_step:end));
    
%     legend('EKF','GPS')
%     legend('EKF','GPS+EKF')
%     legend('UKF','GPS+UKF')
%     legend('EKF','UKF','GPS')
%     legend('GPS+EKF','GPS+UKF')
    legend('GPS+EKF+UWB','GPS+UKF+UWB')
%     legend('GPS','GPS opt')
%     legend('GPS_UKF','GPS')
%     legend('UKF','GPS+UKF')
%     legend('GPS')
%     legend('Chi','GPS')
%     legend('GPS+UKF','GPS+UKF 8 ag')
end

end





