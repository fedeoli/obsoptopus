%OBSERVER SETUP PARAMETERS FILE

global ObserverTest Agent  fh_c r_c h_c  


%ObserverTest.EstimatedPositionInClosedLoop = 1; %1: to use the estimated position in closed loop, otherwise 0
ObserverTest.MaxSteps = min(30,length(time)); %reduced time to this value for shorter simulations
ObserverTest.ObserverON = 1; %Use (1) or not (0) the observer 
ObserverTest.CorrectionON = 1; % UKF: correction step 
ObserverTest.ControlOn = 1; %leave 1, control in the estimated satellite dynamics (to simulate error in the control action, ever)
ObserverTest.ZeroErrors = 1; %to proceed with no measurement errors
ObserverTest.TrasmittedTruePositions = 0; %mainly for debug, assume that true positions and not the estimated ones are transmitted
ObserverTest.TrueInitialPositions = 1;%mainly for debug, assume that true positions are set initially
ObserverTest.UWBoptimizationOn = 1; %GPS is changed in order to minimize a cost funciont with relative distance errors
ObserverTest.StartSoonUWB = 0; %start simulation immediatly optimizing using UWB but with optimization starting condition the one of the GPS: USE IF ALL GPS ARE SEND!
ObserverTest.UWBDropMessages = 0; %some UWB messages are lost using the GILBERT-ELLIOT MODEL (loss bursts)
ObserverTest.GPSDropMessages = 0; %some GPS readings are lost using the GILBERT-ELLIOT MODEL (loss bursts)
ObserverTest.GPSoffset = 0; %try to estimate the constant offset of the GPS (each estimates its own)
ObserverTest.UWBinitialConditionOptimizationIsGPS = 0; %the initial condition of the GPS optimization via UWB has each time the GPS measure as the starting value
ObserverTest.GPSfilterOn = 0; %filtering the x,y,z coordinates of the GPS with three filters
ObserverTest.AddDerivativeEstimation = 2; %=1 obtained by pseudoderivative,= 2 by GPS
ObserverTest.CentralizedOptimization = 0; %optimize all the variables
ObserverTest.OptimizationsEachSampling = 1 ; %please leave 1, normal, more if mutiple optimization among sampling times are performed...more
ObserverTest.UWB_StepsBackInTime = 0; %if 0 is standard, if not we go consider in the J also ObserverTest.UWB_StepsBackInTime 
ObserverTest.EnableSpeedOptimizationWhenSteppingBack = 0; %when ObserverTest.UWB_StepsBackInTime>0 it is possible to optimize the GPS speed, set zero if  ObserverTest.UWB_StepsBackInTime = 0;
ObserverTest.WindowStepMultiplier = 1; %1 by default (>=1), each step back-in-time to compare the actual estimation with past measure is with respect to multiple of the sampling time
ObserverTest.GPS_LVLH2ECI = 1; %1 if GPS noise is in the LHLH coordinates and then is brought back to ECI
ObserverTest.AttitudeObserverOn = 1; %use the attitude observer, otherwise use the true state 
ObserverTest.EstimatedPositionInClosedLoop = 0; %1: to use the estimated position in closed loop, otherwise 0
ObserverTest.version = 104;
ObserverTest.SaveData = 1; %save data if =1 in the corresponding folder as clear below
ObserverTest.FileName = ['Test_UKF_UWB_V' num2str(ObserverTest.version) '_' datestr(now,'yyyymmddTHHMMSS' ,2020) '.mat'];

ObserverTest.Nagents = N_deputy+1; %total number of agents
ObserverTest.Nuwb = ObserverTest.Nagents*(ObserverTest.Nagents-1)/2; %number of distance measurements
ObserverTest.Ndisturbance = 3; %acting on the traslational dynamics of each agent
ObserverTest.Ngps = 3; %number of measure of single GPS (not IMU and others in this first cas
ObserverTest.GPSGaussianCovariance = [15; 10; 5; 0.05; 0.04; 0.02]*1e-3;% GPD error 1-sigma covariance (ref. mail 20190311 from Nardecchia Luca)
%ObserverTest.GPSGaussianCovariance =  [5; 4; 2; 0.05; 0.04; 0.02]*1e-3% GPD error 1-sigma covariance (ref. mail 20190311 from Nardecchia Luca)
ObserverTest.ErrorAmplitudeGPS = 1*max(ObserverTest.GPSGaussianCovariance); %km
ObserverTest.ErrorAmplitudeUWB = 0.15E-3; %km
ObserverTest.IntialConditionPercentage = [1.003 1.002 0.998 1.01 -1.01 -1.03]; %initial percentage error
ObserverTest.IntialConditionAdditive = [0.001, 0.002, -0.001, 0,0,0]; %initial additive error
if(ObserverTest.AddDerivativeEstimation ~= 0)
    ObserverTest.Nspeed = 3;
else
    ObserverTest.Nspeed = 0;
end
ObserverTest.Pi = eye(6)*1E-2;  %P of the UKF estimator, so on below
ObserverTest.Qi = eye(ObserverTest.Ndisturbance)*1E-2;
ObserverTest.Ri = blkdiag(eye(ObserverTest.Ngps)*(30)*1E-3,eye(ObserverTest.Nspeed)*1E-3);%relativa alla misura del GPS e delle velocit? stimate
%ObserverTest.Att_Pi = eye(4+3+2)*1E-3;  %P of the UKF estimator, so on below
%ObserverTest.Att_Qi = eye(ObserverTest.Ndisturbance)*1E-5;
%ObserverTest.Att_Ri = blkdiag(eye(ObserverTest.Ngps)*1E-5,eye(ObserverTest.Nspeed)*1E-5);%relativa alla misura del GPS e delle velocit? stimate
ObserverTest.position_P_reset_aftersamples = 50; %when time counter is mod(, that)==0 the P is reset. Assign equal to zero if no reset is desired

[a,b] = butter(2,2*pi*0.25,'s');
ObserverTest.Filter_position = c2d(tf(a,b),params.time_step,'Tustin');
[ObserverTest.FilterAposition,ObserverTest.FilterBposition,ObserverTest.FilterCposition,ObserverTest.FilterDposition] = ssdata(ObserverTest.Filter_position);
ObserverTest.PseudoDerivative_c = 1;
ObserverTest.PseudoDerivative_d = 1;
ObserverTest.FilterSteps2Regime2EvaluateDerivatives = ObserverTest.PseudoDerivative_d + 0;
ObserverTest.UWBDropMessagesP = 0.2 %probability of losing an UWB message if the previous was sent
ObserverTest.UWBDropMessagesR = 0.85 %probability of sending an UWB message if the previous was lost
ObserverTest.GPSDropMessagesP = 0.05 %probability of losing GPS data if the previous was sent
ObserverTest.GPSDropMessagesR = 0.95 %probability of getting   message if the previous was lost
ObserverTest.optionsObs = optimoptions(@fminunc,  'MaxIter', 100,'display','off'); % 'Algorithm', 'trust-region','SpecifyObjectiveGradient', true,                   

ObserverTest.ExtractGPS = [1 2 3]; %indices of the state vector corresponding to GPS measurements
for k =1:ObserverTest.Nagents,
    ObserverTest.ExchangeDataArrayIndex((k-1)*3+1:(k-1)*3+3) = ObserverTest.ExtractGPS+(k-1)*6;
    ObserverTest.PositionArrayIndex((k-1)*3+1:(k-1)*3+3) = ObserverTest.ExtractGPS+(k-1)*6;
end
ObserverTest.Npassi = length(time)+1;
ObserverTest.Na = (6+ObserverTest.Ndisturbance); %extended state

%TEST parameters
ObserverTest.N_P = 1; %number of changes of P >=2 for the multiple tests
ObserverTest.N_P_min= ObserverTest.Pi*1;
ObserverTest.N_P_max = ObserverTest.Pi*1;
ObserverTest.N_R = 1; %number of changes of R >= 2 for the multiple tests
ObserverTest.N_R_min = ObserverTest.Ri*1;
ObserverTest.N_R_max = ObserverTest.Ri*1;
ObserverTest.expJUWB = 2; %exponent of the distance (UWB) error in J
ObserverTest.expJGPS = 2; %exponent of the distance (GPS) error in J
ObserverTest.expJSIGMA = 2; %exponent of the distance (SIGMA, previous estimated position) error in J
ObserverTest.expJGPSall = 2; %exponent of the distance (GPS all, mean measured GPS from mean estimated position) error in J
ObserverTest.expJGPSpeed = 2; %exponent of GPSpeed error
ObserverTest.expJGPSpeedSigma = 2; %exponent of GPSpeed error
ObserverTest.MeasureMagnifier = 1E3;%the estimation errors are in Km, so when distance becomes below the meter the quadratic/quartic etc cost function of J behaves badly and the error in km is multiplied by this value...
ObserverTest.N_Wuwb_range = [5] %[600 100 ];%[300 150 80 50 30 20 10 5]; %50
ObserverTest.N_Wgps_range = 1%[15 4]%[1000 400 50];%[]
ObserverTest.N_Wsigma_range = 90 % [10 6 4 2] %[10000  6000 3000]%[40 100 250]% [100 50 10];% [ 30000 20000 10000 6000 3000 1000 500 200] %10000
ObserverTest.Weight_GPSall = 3% 30%
ObserverTest.AprioriMismatch4FastConvergence = 1.6E-3% 1.6E-3;
ObserverTest.Weight_Innovation_Fast = 0.2*ObserverTest.N_Wsigma_range;%When the a-priori distance among agents is higher than AprioriMismatch4FastConvergence, use a smaller sigma value to allow faster convergence
ObserverTest.Weight_Innovation_Fast_Reset = 0.05*ObserverTest.Weight_Innovation_Fast;
ObserverTest.Weight_Innovation_Reset_aftersamples =  50; %as P reset but for the optimization of GPS on the Sigma weight
ObserverTest.AprioriMismatch4Projection = 0.1E-2;
ObserverTest.GPSprojectionOnAprioriSphereRadius = 30E-3; %when ObserverTest.AprioriMismatch4FastConvergence is satisfied, then the GPS is projected onto such sphere 
ObserverTest.ThetaLossGPSdata = 0.0; %loss GPS data are substituted with a convex combination of reconstructed position by past data (*ObserverTest.ThetaLossGPSdata) and a-priori estimation *(1-ObserverTest.ThetaLossGPSdata)
ObserverTest.PositionErrorBelowThisAllow2UseEstimatedSpeed = 0.010; %km,  IF ObserverTest.AddDerivativeEstimation = 3
ObserverTest.Weight_UWB_back = ones(1,ObserverTest.UWB_StepsBackInTime+1) %most left entries with most recent data, most right entries with older data - weighting differently back in time the error, then multiplied at runtime for Wuwb (different tests)
ObserverTest.Weight_GPS_back =  ones(1,ObserverTest.UWB_StepsBackInTime+1) %[1 zeros(1,ObserverTest.UWB_StepsBackInTime) ] %ones(1,ObserverTest.UWB_StepsBackInTime+1) % %weighting differently back in time the error
ObserverTest.Weight_GPSall_back = ones(1,ObserverTest.UWB_StepsBackInTime+1) %[1 zeros(1,ObserverTest.UWB_StepsBackInTime)]%ones(1,ObserverTest.UWB_StepsBackInTime+1) %weighting differently back in time the error
ObserverTest.DeadZoneGPS = ObserverTest.ErrorAmplitudeGPS*0.0;
ObserverTest.DeadZoneAllGPS = ObserverTest.ErrorAmplitudeGPS*0.0;
ObserverTest.DeadZoneUWB = ObserverTest.ErrorAmplitudeUWB*0.0;
ObserverTest.UWBOptimizationNoBeforeThan = 10 +ObserverTest.WindowStepMultiplier*ObserverTest.UWB_StepsBackInTime;
ObserverTest.UWBonGPSerror = 15E-1%UWB off if ||GPS - aprioryXgps||>=ObserverTest.UWBonGPSerror
ObserverTest.N_WsigmaOffset = 0;
ObserverTest.StartIntervalWindowPercentage = 0.6; %start of the window for performance evaluation
ObserverTest.EndIntervalWindowPercentage = 0.95; %end of the window for performance evaluation
ObserverTest.tstart = time(1);
ObserverTest.tend = time(end);

ObserverTest.xHatUKF_0 = zeros(ObserverTest.Nagents,6);
ObserverTest.iner_ECI = zeros(ObserverTest.Nagents,6);
ObserverTest.time_step = time_step;
ObserverTest.current_time = time(1);


%ATTITUDE
ObserverTest.AttitudeZeroErrors = 0;
ObserverTest.attitudeTrueInitialPositions = 0;
ObserverTest.AttitudeQ = blkdiag(1,1,1,10,10,10,1,1,1,0.01,0.01,0.01)*1E-6;%by the paper
ObserverTest.AttitudeR = blkdiag(2.7,2.7,2.7,0.012,0.012,0.012)*1E-3; 
ObserverTest.AttitudeP = blkdiag(1,1,1,1,1,1,0.1,0.1,0.1,1,1,1)*1E-3;  %P of the UKF attitude estimator, [quaternion, omega, bias mag, bias gyro]
ObserverTest.MagnetoBias = 0*3E-3*[1;2;3];
ObserverTest.GyroBias = 0*1E-5*[1;3;-1];

ObserverTest.AttitudeInitialConditionAdditive_Sigma_EulerAngles = 10*pi/180; %initial percentage error on angles
ObserverTest.AttitudeInitialConditionAdditive_Sigma_Omega = 1E-4;
ObserverTest.AttitudeInitialConditionAdditive_Sigma_GyrosBias = 0*1E-5;
ObserverTest.AttitudeInitialConditionAdditive_Sigma_MagnetoBias = 0*1E-1;
ObserverTest.NumberOfAttitudeDifferentInitialConditions = 1;%simulations are performed a number of times with radn initial conditions havinf the sigma just above

ObserverTest.headingMagDisplacement = 0*pi/180;
ObserverTest.tiltMagDisplacement = 0*pi/180;
ObserverTest.MagSigmaElevation = 0*pi/180;%sigma of noise in spherical coordinates for the Elevation component
ObserverTest.MagSigmaAzimut = 0*pi/180;%sigma of noise in spherical coordinates for the Azimutal component
ObserverTest.GyroGaussianCovariance = 3*[1 1 1]*1E-4;
ObserverTest.headingGyroDisplacement = 0*pi/180;
ObserverTest.tiltGyroDisplacement = 0*pi/180;
ObserverTest.EulerAngleNoiseOnMag = 0.01;

ObserverTest.attitudeBackSamplesWindow = 10; %how many samples to fit (window size for data matching)
ObserverTest.attitudeBackIntersample = 3;%number of steps among samples
ObserverTest.attitudeWnormquaternion = 1E5;
ObserverTest.attitude_Km = 1;
ObserverTest.attitude_Kmbis = 1;
ObserverTest.attitude_Ka = 0*3;
ObserverTest.attitude_Kp = 1;
ObserverTest.attitude_Ki = 0*0.3;
ObserverTest.attitude_Komega = 2;

%ObserverTest.SwitchOnAttitudeObserverAfterSteps = 10;

% Data for iterative test are created
ObserverTest.TotalRuns = ObserverTest.N_P*ObserverTest.N_R*length(ObserverTest.N_Wuwb_range)*length(ObserverTest.N_Wgps_range)*length(ObserverTest.N_Wsigma_range)*ObserverTest.NumberOfAttitudeDifferentInitialConditions;
ObserverTest.MaxError = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);
ObserverTest.Mean = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);
ObserverTest.Sigma = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);
ObserverTest.AttitudeMaxError = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);
ObserverTest.AttitudeMean = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);
ObserverTest.AttitudeSigma = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);

ObserverTest.TransientTime = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);
ObserverTest.UnfeasableR = zeros(ObserverTest.TotalRuns,ObserverTest.Nagents);


kk = 0;
Pindexes = zeros(1,ObserverTest.TotalRuns);
Rindexes = zeros(1,ObserverTest.TotalRuns);
for np=1:ObserverTest.N_P,
    for nr=1:ObserverTest.N_R,
        for nuwb=1:length(ObserverTest.N_Wuwb_range),
            for ngps=1:length(ObserverTest.N_Wgps_range),
                for nsigma=1:length(ObserverTest.N_Wsigma_range),
                    for aic=1: ObserverTest.NumberOfAttitudeDifferentInitialConditions, 
                        
                    kk = kk+1;
                    ObserverTest.AttitudeInitialConditionAdditive_EulerAngles = randn(1,3)*ObserverTest.AttitudeInitialConditionAdditive_Sigma_EulerAngles; %initial percentage error on angles
                    ObserverTest.AttitudeInitialConditionAdditive_Omega = randn(1,3)*ObserverTest.AttitudeInitialConditionAdditive_Sigma_Omega;
                    ObserverTest.AttitudeInitialConditionAdditive_Bias = [randn(1,3)*ObserverTest.AttitudeInitialConditionAdditive_Sigma_MagnetoBias,randn(1,3)*ObserverTest.AttitudeInitialConditionAdditive_Sigma_GyrosBias];
    
                    theta = (ObserverTest.N_P-(np-1))/ObserverTest.N_P;
                    Pindexes(np) = theta;
                    Pi_temp = ObserverTest.N_P_max*theta+ (1- theta)*ObserverTest.N_P_min; %convex combination
                    [n1,n2] = size(Pi_temp);
                    ObserverTest.PiAll(kk,:) = reshape(Pi_temp,[1,n1*n2]);
                    theta = (ObserverTest.N_R-(nr-1))/ObserverTest.N_R;
                    Rindexes(nr) = theta;
                    Ri_temp = ObserverTest.N_R_max*theta + (1- theta)*ObserverTest.N_R_min;
                    [n1,n2] = size(Ri_temp);
                    ObserverTest.RiAll(kk,:) = reshape(Ri_temp,[1,n1*n2]);
                    ObserverTest.Gains(:,kk) = [ObserverTest.N_Wuwb_range(nuwb);ObserverTest.N_Wgps_range(ngps);ObserverTest.N_Wsigma_range(nsigma);];
                    end
                end
            end
        end
    end
end














