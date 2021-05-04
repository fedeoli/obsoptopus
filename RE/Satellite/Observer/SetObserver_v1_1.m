function [ObserverTest, Agent] = SetObserver_v1_1(iteration,lengthtime,satellites_iner_ECI,satellites_attitude,params)
%%load all observer parameters
% global ObserverTest Agent

S = load('ObserverTest.mat');
ObserverTest = S.ObserverTest;

%%%%%%%% 05/10 %%%%%%%
ObserverTest.Nagents = params.Ndeputy + 1;

ObserverTest.ExchangeDataArrayIndex = zeros(1,3*ObserverTest.Nagents);
for i = 1:ObserverTest.Nagents 
   ObserverTest.ExchangeDataArrayIndex(1 + 3*(i-1):3 + 3*(i-1)) = [1 2 3] + 6*(i-1);
end

ObserverTest.PositionDataArrayIndex = zeros(1,3*ObserverTest.Nagents);
ObserverTest.PositionArrayIndex = zeros(1,3*ObserverTest.Nagents);
for i = 1:ObserverTest.Nagents 
   ObserverTest.PositionArrayIndex(1 + 3*(i-1):3 + 3*(i-1)) = [1 2 3] + 6*(i-1);
end

%%%%%% FLAG + NOISE SETUP %%%%%%%
% Position Observer flags
ObserverTest.ObserverON = 1;
% GPS flags
ObserverTest.GPSopt_flag = 0;
ObserverTest.StartSoonUWB = 0; % set to 1 if only GPS optimization is set
ObserverTest.UWBOptimizationNoBeforeThan = 60; % used when KF is enabled
% KF flags
ObserverTest.KF_flag = 1;
ObserverTest.KF_pos = 'UKF';
% Attitude observer flags
ObserverTest.AttitudeObserverOn = 0;

% GPS optimization parameters - geometric method
ObserverTest.theta = 0.03;
ObserverTest.beta = 0.5;
ObserverTest.check_distance = 0;
ObserverTest.projection = 'Chi';

% Measurement bias
ObserverTest.MagnetoBias = 5e-7*randn(3,1) + 1e-6;
ObserverTest.GyroBias = 5e-2*randn(3,1) + 1e-1;
% ObserverTest.ErrorAmplitudeGPS = 5e-3;    PER ORA INUTILE (VEDI OBSERVER_MEASUREMENTS_V1_3)
ObserverTest.GPSGaussianCovariance = [15; 10; 5; 0.05; 0.04; 0.02]*1e-3;
ObserverTest.ErrorAmplitudeUWB = 5e-4;

%%%%% COVARIANCE ATTITUDE %%%%%
% ObserverTest.AttitudeQ = diag([1e-6 1e-6 1e-6 1e-5 1e-5 1e-5 1e-6 1e-6 1e-6 1e-8 1e-8 1e-8 1e-6 1e-6 1e-6]);
% ObserverTest.AttitudeP = diag([1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-4 1e-4 1e-4 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3]);
ObserverTest.AttitudeQ = 1e-6*eye(3*ObserverTest.Nagents);
ObserverTest.AttitudeP = 1e-3*eye(3*ObserverTest.Nagents);

%%%%% COVARIANCE POSITION %%%%%
ObserverTest.Pi = eye(6)*1E-2;  %P of the UKF estimator, so on below
ObserverTest.Qi = eye(ObserverTest.Ndisturbance)*1E-2;
ObserverTest.Ri = blkdiag(eye(ObserverTest.Ngps)*(30)*1E-3,eye(ObserverTest.Nspeed)*1E-3);
%%%%%%%%%%%%%%%%%%%%%%

% Initial condition noise
ObserverTest.IntialConditionAdditive = 1*ObserverTest.IntialConditionAdditive;
ObserverTest.IntialConditionPercentage = 1*ObserverTest.IntialConditionPercentage;

% Packet loss
ObserverTest.UWBDropMessages = 0;
ObserverTest.UWBDropMessagesP = 0.3; %probability of losing an UWB message if the previous was sent
ObserverTest.UWBDropMessagesR = 0.85; %probability of sending an UWB message if the previous was lost

ObserverTest.GPSDropMessages = 0;
ObserverTest.ThetaLossGPSdata = 0;
ObserverTest.GPSDropMessagesP = 0.3; %probability of losing GPS data if the previous was sent
ObserverTest.GPSDropMessagesR = 0.95; %probability of getting   message if the previous was lost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObserverTest.MaxSteps = lengthtime;
ObserverTest.Npassi = lengthtime;
ObserverTest.Ntest = iteration;
[a,b] = butter(2,2*pi*0.25,'s');
ObserverTest.Filter_position = c2d(tf(a,b),params.time_step,'Tustin');
[ObserverTest.FilterAposition,ObserverTest.FilterBposition,ObserverTest.FilterCposition,ObserverTest.FilterDposition] = ssdata(ObserverTest.Filter_position);
ObserverTest.Npassi = lengthtime;
ObserverTest.Na = (6+ObserverTest.Ndisturbance); %extended state
ObserverTest.time_step = params.time_step;


%Settings variables before the estimation process selecting the opportune
%initial estimates, the same for all runs
InitialConditionPercentageMatrix = diag(ObserverTest.IntialConditionPercentage);
ObserverTest.satellites_iner_ECI_0 = satellites_iner_ECI;
ObserverTest.satellites_attitude_0 = satellites_attitude;
for k = 1:ObserverTest.Nagents
    ObserverTest.iner_ECI_0(k,:) = ObserverTest.satellites_iner_ECI_0(1+(k-1)*6:6+(k-1)*6); %initial position state
    ObserverTest.attitude_0(k,:) = ObserverTest.satellites_attitude_0(1+(k-1)*7:7+(k-1)*7); %initial attitude
    ObserverTest.attitude_0(k,1:4) = quatnormalize(ObserverTest.attitude_0(k,1:4));
    ObserverTest.attitude_xHatUKF_0(k,:) = zeros(4+3+6,1);
    if(ObserverTest.TrueInitialPositions==1)
        ObserverTest.xHatUKF_0(k,:) = ObserverTest.iner_ECI_0(k,:);
    else
        ObserverTest.xHatUKF_0(k,:) = InitialConditionPercentageMatrix*ObserverTest.iner_ECI_0(k,:)'+(ObserverTest.IntialConditionAdditive*2*(rand-0.5))'; %estimated initial state
    end
    if(ObserverTest.attitudeTrueInitialPositions==1)
        ObserverTest.attitude_xHatUKF_0(k,1:7) = ObserverTest.attitude_0(k,:);
        ObserverTest.attitude_xHatUKF_0(k,8:end) = [ObserverTest.MagnetoBias;ObserverTest.GyroBias]; 
    else
        R = quat2dcm(ObserverTest.attitude_0(k,1:4));
        R = angle2dcm(ObserverTest.AttitudeInitialConditionAdditive_EulerAngles(1),ObserverTest.AttitudeInitialConditionAdditive_EulerAngles(2),ObserverTest.AttitudeInitialConditionAdditive_EulerAngles(3))*R;
        ObserverTest.attitude_xHatUKF_0(k,1:4) = dcm2quat(R);
        ObserverTest.attitude_xHatUKF_0(k,5:7) = (ObserverTest.AttitudeInitialConditionAdditive_Omega + ObserverTest.attitude_0(k,5:7))';
        ObserverTest.attitude_xHatUKF_0(k,8:end) = (ObserverTest.AttitudeInitialConditionAdditive_Bias + [ObserverTest.MagnetoBias;ObserverTest.GyroBias]'); 
    end
end

%Switching off the warnings risen by the optimization procedure
warning('off','optim:fminunc:SwitchingMethod');

ObserverTest.AdjacencyMatrix = zeros(ObserverTest.Nagents,ObserverTest.Nagents);
ObserverTest.satellites_iner_ECI_alltimes = zeros(6*ObserverTest.Nagents,ObserverTest.Npassi);
ObserverTest.estimated_satellites_iner_ECI = zeros(6*ObserverTest.Nagents,ObserverTest.Npassi);
ObserverTest.satellites_iner_ECI_alltimes(:,1) = ObserverTest.satellites_iner_ECI_0 ;
ObserverTest.satellites_iner_ECI = ObserverTest.satellites_iner_ECI_0 ;
ObserverTest.params_u_alltimes = zeros(3*(ObserverTest.Nagents-1),ObserverTest.Npassi); 
ObserverTest.estimated_satellites_attitude = zeros(7*ObserverTest.Nagents,ObserverTest.Npassi);

ObserverTest.Weight_UWB = ObserverTest.Gains(1,ObserverTest.Ntest);
ObserverTest.Weight_GPS = ObserverTest.Gains(2,ObserverTest.Ntest);
ObserverTest.Weight_Innovation = ObserverTest.Gains(3,ObserverTest.Ntest);
ObserverTest.Weight_GPSall = ObserverTest.Weight_GPSall;

ObserverTest.GPS_alltimes = zeros(3*ObserverTest.Nagents,ObserverTest.Npassi);
ObserverTest.GPS_optimized_alltimes = zeros(3*ObserverTest.Nagents,ObserverTest.Npassi);

ObserverTest.Weight_UWB_back = ObserverTest.Weight_UWB_back*ObserverTest.Weight_UWB;
ObserverTest.Weight_GPS_back =  ObserverTest.Weight_GPS_back*ObserverTest.Weight_GPS;
ObserverTest.Weight_GPSall_back = ObserverTest.Weight_GPSall_back*ObserverTest.Weight_GPSall;

% define for C coder
Agent = struct();


for k = 1:ObserverTest.Nagents
    
    Agent(k).iner_ECI = zeros(6,ObserverTest.Npassi); %true state of each satellite
    Agent(k).xHatUKF = zeros(6,ObserverTest.Npassi); %Estimated state
    Agent(k).iner_ECI_EstimationError = zeros(6,ObserverTest.Npassi); %estimation errors
    Agent(k).NormEstimationErrorXYZ = zeros(ObserverTest.Npassi); %estimation error norm for x,y,z
    Agent(k).GotData = zeros((ObserverTest.UWB_StepsBackInTime+1)*ObserverTest.WindowStepMultiplier,length(ObserverTest.ExchangeDataArrayIndex)); %data exchanded in the constellation
    Agent(k).AgentsEstimatedPositions = Agent(k).GotData(1,:); %self-estimated positions of the other agents
    Agent(k).iner_ECI(:,1) = ObserverTest.iner_ECI_0(k,:); %initial state
    Agent(k).xHatUKF(:,1) = ObserverTest.xHatUKF_0(k,:);
    Agent(k).xaUKF = [Agent(k).xHatUKF(:,1); zeros(ObserverTest.Ndisturbance,1)]; %Extended state (single step)
    Agent(k).iner_ECI_EstimationError(:,1) = Agent(k).iner_ECI(:,1) - Agent(k).xHatUKF(:,1);
    Agent(k).NormEstimationErrorXYZ(1) = norm(Agent(k).iner_ECI_EstimationError(1:3,1));
    Agent(k).xaUKF(1:6,1) = Agent(k).xHatUKF(:,1);    % the estimate of r is subtracted by params.Re
    Agent(k).P = reshape(ObserverTest.PiAll(ObserverTest.Ntest,:),size(ObserverTest.Pi));
    Agent(k).Q = ObserverTest.Qi;
    Agent(k).R = reshape(ObserverTest.RiAll(ObserverTest.Ntest,:),size(ObserverTest.Ri));
    Agent(k).UWBcorrectionActive = zeros(1,ObserverTest.Npassi);
    ObserverTest.estimated_satellites_iner_ECI(1+(k-1)*6:6+(k-1)*6,1) = Agent(k).xHatUKF(:,1);
    Agent(k).GPS = zeros(3,ObserverTest.Npassi);
    Agent(k).GPSpeed = zeros(3,ObserverTest.Npassi);
    Agent(k).GPSopt = zeros(3,ObserverTest.Npassi);
    Agent(k).GPSfilter = zeros(3*length(ObserverTest.FilterBposition),ObserverTest.Npassi); %state of the filter for the GPS along x,y,z
    Agent(k).GPSfiltered = zeros(3,ObserverTest.Npassi); %state of the filter for the GPS along x,y,z
    Agent(k).GPSderivative = zeros(3,ObserverTest.Npassi); %state of the filter for the GPS along x,y,z
    Agent(k).GPSpeedOpt = Agent(k).GPSpeed;
    Agent(k).GPSoffset = [k/4; 0; -k/4]; %GPS offset inizialization
    Agent(k).GPSoffsetEstimated = zeros(3,ObserverTest.Npassi); %GPS offset estimated
    if(ObserverTest.ZeroErrors==1)
        myGPS = Agent(k).iner_ECI(ObserverTest.ExtractGPS,1);
        myGPSpeed = Agent(k).iner_ECI(4:6,1);
    else
        noise = ObserverTest.GPSGaussianCovariance.*randn(6,1);
        myGPS = noise(1:3)+ Agent(k).iner_ECI(ObserverTest.ExtractGPS,1);
        myGPSpeed = noise(4:6) + Agent(k).iner_ECI(4:6,1);
        %myGPS = (ErrorAmplitudeGPS*2*(0.5*ones(1,Ngps)-rand(1,Ngps)))'...
        %   + Agent(k).iner_ECI(ExtractGPS,1); %GPS single satellite measurement: (x,y,z)
    end
    Agent(k).GPS(:,1) = myGPS;
    Agent(k).GPSpeed(:,1) = myGPSpeed;
    Agent(k).GPSopt(:,1) = myGPS;
    Agent(k).GPSpeedOpt(:,1) = Agent(k).GPSpeed(:,1);
    for jj=1:3,
        Agent(k).GPSfilter([1 2] + (jj-1)*2,1) = pinv(eye(2)-ObserverTest.FilterAposition)*ObserverTest.FilterBposition*myGPS(jj);
    end
    Agent(k).SuccessfullyReceivedData = ones((ObserverTest.UWB_StepsBackInTime+1)*ObserverTest.WindowStepMultiplier,ObserverTest.Nagents);
    Agent(k).SuccessfullyReadGPS = ones(1,(ObserverTest.UWB_StepsBackInTime+1)*ObserverTest.WindowStepMultiplier);
    Agent(k).MeasuredDistances = zeros(ObserverTest.Npassi,ObserverTest.Nagents);
    Agent(k).MismatchAprioriDistances = zeros(1,ObserverTest.Npassi);
    Agent(k).MismatchAposterioriDistances = zeros(1,ObserverTest.Npassi);
    Agent(k).ErrorAposterioriDistances = zeros(1,ObserverTest.Npassi);
    Agent(k).Pmeno = ObserverTest.Pi*0;
    Agent(k).xHatMeno = zeros(6,1);
    Agent(k).xHatUKFmultiple = Agent(k).xHatUKF(:,1);
    
    
    Agent(k).attitude = zeros(7,ObserverTest.Npassi);
    Agent(k).attitude_xHatUKF = zeros(7+6,ObserverTest.Npassi);
    Agent(k).magneto = zeros(3,ObserverTest.Npassi); %measurements from magnetometer
    Agent(k).magneto2 = zeros(3,ObserverTest.Npassi); %measurements from magnetometer
    Agent(k).gyro = zeros(3,ObserverTest.Npassi); %measurements from gyro
    Agent(k).magnetoBias = zeros(3,ObserverTest.Npassi); 
    Agent(k).gyroBias = zeros(3,ObserverTest.Npassi);
    Agent(k).attitude_error = zeros(7,ObserverTest.Npassi);
    Agent(k).magyrobias_error = zeros(6,ObserverTest.Npassi);
    
    Agent(k).magnetoBias(:,1) = ObserverTest.MagnetoBias;
    Agent(k).gyroBias(:,1) = ObserverTest.GyroBias;
    Agent(k).attitude(:,1) = ObserverTest.attitude_0(k,:)';
    [MagTemp,MagTemp2,GyrosTemp] = AttitudeObserver_GetMeasures_v2_3(Agent(k).iner_ECI(1:3,1)',Agent(k).attitude(1:4,1)',...
                                                                                                        Agent(k).attitude(5:7,1)',Agent(k).magnetoBias(:,1)',Agent(k).gyroBias(:,1)',ObserverTest.RPYbetweenMagSensors,0,ObserverTest);
    Agent(k).magneto(:,1) = MagTemp;
    Agent(k).magneto2(:,1) = MagTemp2;
    Agent(k).gyros(:,1) = GyrosTemp;
    Agent(k).attitude_xHatUKF(:,1) = ObserverTest.attitude_xHatUKF_0(k,:)';
    
    %%%%% set R_ECI2Body - 06/10 %%%%%
    Agent(k).R_ECI2Body = store_dcm(Agent(k).attitude_xHatUKF(:,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Agent(k).attitudeP = ObserverTest.AttitudeP;
    Agent(k).attitude_error(1:4,1) = quatnormalize(quatmultiply(quatinv(Agent(k).attitude(1:4,1)'), Agent(k).attitude_xHatUKF(1:4,1)')); %by CHALLA (8) but we can also swap the inverse among vectors in the quatmultiply to have the standard true-estimate
    Agent(k).attitude_error(5:7,1) = Agent(k).attitude(5:7,1) - Agent(k).attitude_xHatUKF(5:7,1) ;
    Agent(k).magyrobias_error(:,1) = [ObserverTest.MagnetoBias;ObserverTest.GyroBias] - Agent(k).attitude_xHatUKF(8:end,1);
    
    Agent(k).stima(:,1) = zeros(3*3,1); 
    Agent(k).stima(1:3,1) =  Agent(k).magneto(:,1);
    Agent(k).stima(4:6,1) =Agent(k).magneto2(:,1);
    
    
end

%initialization of the data collected by the network
for k = 1:ObserverTest.Nagents
    for jj=1:ObserverTest.Nagents
        if(ObserverTest.TrasmittedTruePositions==1)
            temp = Agent(jj).iner_ECI(1:3,1);
        else
            temp = Agent(jj).xHatUKF(1:3,1);
        end
        Agent(k).GotData(1,(jj-1)*3+1:(jj-1)*3+3) =  temp;
        Agent(k).AgentsEstimatedPositions((jj-1)*3+1:(jj-1)*3+3) = temp;
    end
end
ObserverTest.CSI = zeros(ObserverTest.Na,2*ObserverTest.Na+1); % i sigma points
ObserverTest.CSI_X_meno = zeros(6,2*ObserverTest.Na+1);
ObserverTest.W = zeros(2*ObserverTest.Na+1,1); % i loro pesi
ObserverTest.UestimatorStory = zeros(3*(ObserverTest.Nagents-1),ObserverTest.Npassi);

% define for C code
coe = zeros(1,6);

% assignment
coe(1:6) = rv2coe_V1_1(Agent(1).xHatUKF(1:3,1),Agent(1).xHatUKF(4:6,1),params.mi);
ObserverTest.EstimatedChiefCoe = coe;

if params.Nagents ~= 1
    ObserverTest.DVsizeMemo = size(params.DV(:, :, 1)); 
else
    ObserverTest.DVsizeMemo = 0; 
end


