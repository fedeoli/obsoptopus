%% observer - storage
global ObserverTest

i = length(time);

%Store true and estimated data for all experiments (test)
ObserverTest.AllSimulation(ObserverTest.Ntest).iner_ECI = ObserverTest.satellites_iner_ECI_alltimes;
ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_iner_ECI = ObserverTest.estimated_satellites_iner_ECI;
ObserverTest.AllSimulation(ObserverTest.Ntest).attitude = satellites_attitude_out;
ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_attitude = ObserverTest.estimated_satellites_attitude;

% start_step = 1;
start_step = max(1,floor(ObserverTest.StartIntervalWindowPercentage*(i)));

% end_step = ObserverTest.Npassi;
end_step = floor(ObserverTest.EndIntervalWindowPercentage*(i));

window_interval = start_step:1:end_step;
ObserverTest.window_interval = window_interval;
ObserverTest.window(ObserverTest.Ntest,:) = time([window_interval(1),window_interval(end)]);

for k=1:ObserverTest.Nagents
    ObserverTest.Mean_truepos(ObserverTest.Ntest,k) = norm(mean(Agent(k).iner_ECI(1:3,window_interval),2)); 
    ObserverTest.mean_measure(ObserverTest.Ntest,k) = norm(mean(Agent(k).xHatUKF(1:3,window_interval),2));
    ObserverTest.error_ratio(ObserverTest.Ntest,k) = ObserverTest.mean_measure(ObserverTest.Ntest,k)/ObserverTest.Mean_truepos(ObserverTest.Ntest,k);
    
    ObserverTest.Mean(ObserverTest.Ntest,k) = mean(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.Mean_rel(ObserverTest.Ntest,k) = mean(Agent(k).NormEstimationErrorXYZ(window_interval)/ObserverTest.Mean_truepos(ObserverTest.Ntest,k))*100;
    
    ObserverTest.min_err(ObserverTest.Ntest,k) = min(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.min_rel_err(ObserverTest.Ntest,k) = min(Agent(k).NormEstimationErrorXYZ(window_interval)/ObserverTest.Mean_truepos(ObserverTest.Ntest,k))*100;
    
    ObserverTest.Sigma(ObserverTest.Ntest,k) =  std(Agent(k).NormEstimationErrorXYZ(window_interval));    
    ObserverTest.Sigma_rel(ObserverTest.Ntest,k) = std(Agent(k).NormEstimationErrorXYZ(window_interval)/ObserverTest.Mean_truepos(ObserverTest.Ntest,k))*100;
    
    ObserverTest.maxError(ObserverTest.Ntest,k) =  max(Agent(k).NormEstimationErrorXYZ(window_interval));
%     ObserverTest.TransientTime(ObserverTest.Ntest,k) = find(Agent(k).NormEstimationErrorXYZ <= ObserverTest.Mean(ObserverTest.Ntest,k),1);
%     ObserverTest.MeanDistances(ObserverTest.Ntest,k) = mean(Agent(k).ErrorAposterioriDistances(window_interval));
%     ObserverTest.SigmaDistances(ObserverTest.Ntest,k) =  std(Agent(k).ErrorAposterioriDistances(window_interval));

% NB: not considering Ntest (assuming doing only one)
    ObserverTest.Mean_sign(k,:) = mean(Agent(k).iner_ECI_EstimationError(1:3,window_interval),2);
    ObserverTest.Sigma_sign(k,:) =  std(Agent(k).iner_ECI_EstimationError(1:3,window_interval),0,2);
    
    ObserverTest.Gpserror(k,:,:) = Agent(k).iner_ECI(1:3,:) - Agent(k).GPS;
    ObserverTest.Mean_Gpserror(k,:) = mean(ObserverTest.Gpserror(k,:,window_interval),3);
    ObserverTest.Sigma_Gpserror(k,:) = std(ObserverTest.Gpserror(k,:,window_interval),0,3);
end

% index creation
% voglio far vedere che la sigma cala con la stima
temp_sigma = ObserverTest.Sigma_sign./ObserverTest.Sigma_Gpserror;
% voglio far vedere che le due media sono circa uguali, sempre attorno a 0
temp_mean = abs(ObserverTest.Mean_sign - ObserverTest.Mean_Gpserror);
for k = 1:ObserverTest.Nagents
   ObserverTest.sigma_gain(1,k) = mean(temp_sigma(k,:),2); 
   ObserverTest.mean_gain(1,k) = mean(temp_mean(k,:),2);
end

if(ObserverTest.SaveData)
    if(ObserverTest.CentralizedOptimization==1)
        save(['Observer/CentralizedTest/' ObserverTest.FileName],'ObserverTest');
    else
        save(['Observer/DecentralizedTest/' ObserverTest.FileName],'ObserverTest');
    end
end

for k = 1:ObserverTest.Nagents
    Agent(k).q_Euler = zeros(ObserverTest.Npassi,3);
    Agent(k).q_est_Euler = zeros(ObserverTest.Npassi,3);
     Agent(k).delq_Euler = zeros(ObserverTest.Npassi,3);
    for z = 1:ObserverTest.Npassi
        Agent(k).q_Euler(z,:) = quat2eul(Agent(k).attitude(1:4,z)', 'ZYX')'; 
        Agent(k).q_est_Euler(z,:) = quat2eul(Agent(k).attitude_xHatUKF(1:4,z)', 'ZYX')';
        Agent(k).delq_Euler(z,:) = Agent(k).q_Euler(z,:) - Agent(k).q_est_Euler(z,:);
    end
end

for i = 1:ObserverTest.Nagents
    Chi = ObserverTest.AllSimulation(ObserverTest.Ntest).iner_ECI(1+6*(i-1):3+6*(i-1),:);
    
    % GPS trajectory
    Chi_est_GPS = Agent(i).GPSopt;
    
    % GPS+ KF trajectory
    Chi_est_KF = Agent(i).xHatUKF;
    
    % compute the norm
    for z = 1:ObserverTest.Npassi
       Chi_norm(z) = norm(Chi(:,z));
       Chi_est_GPS_norm(z) = norm(Chi_est_GPS(:,z));
       Chi_est_KF_norm(z) = norm(Chi_est_KF(:,z));
    end
    
    % GPS and KF error
    traj_error_GPS(:,i) = abs(Chi_norm-Chi_est_GPS_norm);
    traj_error_KF(:,i) = abs(Chi_norm-Chi_est_KF_norm);
end

% error mean on attitude
ObserverTest.MeanAttitude = zeros(ObserverTest.Nagents,3);
ObserverTest.SigmaAttitude = zeros(ObserverTest.Nagents,3);
for k = 1:ObserverTest.Nagents
    ObserverTest.MeanAttitude(k,:) = mean(Agent(k).delq_Euler(window_interval,:),1)*180/pi;
    ObserverTest.SigmaAttitude(k,:) = std(Agent(k).delq_Euler(window_interval,:),1)*180/pi;
end


% get sum of error + mean and sigma
traj_error_GPS_sum = sum(traj_error_GPS,2);
mean_traj_error_GPS_sum = mean(traj_error_GPS,2);
sigma_traj_error_GPS_sum = var(traj_error_GPS,0,2);

traj_error_KF_sum = sum(traj_error_KF,2);
mean_traj_error_KF_sum = mean(traj_error_KF,2);
sigma_traj_error_KF_sum = var(traj_error_KF,0,2);