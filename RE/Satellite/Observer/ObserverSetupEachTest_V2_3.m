%Settings variables before the estimation process selecting the opportune
%gains 


%initial estimates, the same for all runs
InitialConditionPercentageMatrix = diag(ObserverTest.IntialConditionPercentage);
ObserverTest.satellites_iner_ECI_0 = satellites_iner_ECI;
ObserverTest.satellites_attitude_0 = satellites_attitude;
for k = 1:ObserverTest.Nagents,
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
        [r,p,y] = quat2angle(ObserverTest.attitude_0(k,1:4));
        temp = ObserverTest.AttitudeInitialConditionAdditive_EulerAngles + [r,p,y];
        ObserverTest.attitude_xHatUKF_0(k,1:4) = angle2quat(temp(1),temp(2),temp(3));
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




for k = 1:ObserverTest.Nagents,
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
    Agent(k).gyro = zeros(3,ObserverTest.Npassi); %measurements from gyro
    Agent(k).magnetoBias = zeros(3,ObserverTest.Npassi); 
    Agent(k).gyroBias = zeros(3,ObserverTest.Npassi);
    Agent(k).attitude_error = zeros(7,ObserverTest.Npassi);
    Agent(k).bias_error = zeros(6,ObserverTest.Npassi);
    
    Agent(k).magnetoBias(:,1) = ObserverTest.MagnetoBias;
    Agent(k).gyroBias(:,1) = ObserverTest.GyroBias;
    Agent(k).attitude(:,1) = ObserverTest.attitude_0(k,:)';
    [MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(Agent(k).iner_ECI(1:3,1)',Agent(k).attitude(1:4,1)',...
                                                                                                        Agent(k).attitude(5:7,1)',Agent(k).magnetoBias(:,1)',Agent(k).gyroBias(:,1)',0);
    Agent(k).magneto(:,1) = MagTemp;
    Agent(k).gyros(:,1) = GyrosTemp;
    Agent(k).attitude_xHatUKF(:,1) = ObserverTest.attitude_xHatUKF_0(k,:)';   
    
    Agent(k).attitudeP = ObserverTest.AttitudeP;
    Agent(k).attitude_error(1:4,1) = quatnormalize(quatmultiply(quatinv(Agent(k).attitude(1:4,1)'), Agent(k).attitude_xHatUKF(1:4,1)')); %by CHALLA (8) but we can also swap the inverse among vectors in the quatmultiply to have the standard true-estimate
    Agent(k).attitude_error(5:7,1) = Agent(k).attitude(5:7,1) - Agent(k).attitude_xHatUKF(5:7,1) ;
    Agent(k).magyrobias_error(:,1) = [ObserverTest.MagnetoBias;ObserverTest.GyroBias] - Agent(k).attitude_xHatUKF(8:end,1);
                
end

%initialization of the data collected by the network
for k = 1:ObserverTest.Nagents,
    for jj=1:ObserverTest.Nagents,
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
ObserverTest.StateChief4Control = zeros(6,1); %it will change size...

