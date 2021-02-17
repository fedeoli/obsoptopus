function Observer_Measurements_V1_3(satellites_iner_ECI, deputy_rel_LVLH, params,i)

%INPUT: 
% satellites_iner: Array (6*N x 1), with N equal to the number of satellites. It is built in the following way: [pos1; vel1; pos2; vel2; pos3; vel3;
%                    .... ], where the position and velocities are must be expressed in ECI reference frame. The first six components (pos1; vel1) are 
%                    the inertial coordinates of the chief satellite. 
% deputy_rel_LVLH: deputy rLVLH coordinates with respect to chief
% params:  object defined by La Sapienza containing general costants,  simulation parameters, control torques, desired references...
% i: previous time step index (the estimates are computed at time t = (i+1)*time_step), IT IS ASSUMED THAT time_step_att = time_step  (position).

%OTUPUT: 
%The global (objects) variables Agent and ObserverTest are updated (fields associated to measurements)


global  fh_c r_c h_c Agent ObserverTest

% Set the AdjacencyMatrix (distence matrices, simmetric (zeros the lower triangular part, non zero the upper triangular one))
%the distance measurament and data exchanged is assumed to be "instantaneous"
ObserverTest.AdjacencyMatrix = setAdjacencyMatrixNorm(satellites_iner_ECI(ObserverTest.PositionArrayIndex),ObserverTest.Nagents);

%AdjacencyMatrix obtained by the apriori sent GPS
ObserverTest.AprioriAdjacencyMatrix =  setAdjacencyMatrixNorm(ObserverTest.APrioriEstimationXYZ,ObserverTest.Nagents);
                        

%relative distances
ObserverTest.MeasuredDistances = zeros(ObserverTest.Nagents);
for n = 1:ObserverTest.Nagents,
    %Measured distances: adding the measurement error to UWB (relative distances)
    for jj = n+1:ObserverTest.Nagents,
        if(ObserverTest.ZeroErrors == 1)
            ObserverTest.MeasuredDistances(n,jj) = ObserverTest.AdjacencyMatrix(n,jj);
        else
            ObserverTest.MeasuredDistances(n,jj) = ObserverTest.ErrorAmplitudeUWB*2*(0.5*(1-rand)) + ObserverTest.AdjacencyMatrix(n,jj);
        end
    end
end
ObserverTest.MeasuredDistances = ObserverTest.MeasuredDistances + ObserverTest.MeasuredDistances';

for n = 1:ObserverTest.Nagents,
    Agent(n).MeasuredDistances(i+1,:) = ObserverTest.MeasuredDistances(n,:);
    %Some GPS readings might fail, keep track
    Agent(n).SuccessfullyReadGPS(1:end-1) = Agent(n).SuccessfullyReadGPS(2:end);    
    
    if(ObserverTest.GPSDropMessages == 1)% &&  i+1> (ObserverTest.UWB_StepsBackInTime+1)*ObserverTest.WindowStepMultiplier)
        Agent(n).SuccessfullyReadGPS(end) = GPSFineReading_v1(Agent(n).SuccessfullyReadGPS(end),ObserverTest.GPSDropMessagesP,ObserverTest.GPSDropMessagesR);
    else
        Agent(n).SuccessfullyReadGPS(end) = 1;
    end
end

%GNSS measurements
satellite_pos_ECI(:,1) = satellites_iner_ECI(1:3);                                                                                                          % Satellite's inertial position vector expressed in ECI (chief)
satellite_vel_ECI(:,1) = satellites_iner_ECI(4:6);                                                                                                          % Satellite's inertial velocity vector expressed in ECI (chief)
coe = rv2coe_V1_1(satellite_pos_ECI(:,1), satellite_vel_ECI(:,1), params.mi);
chief_omega_LVLH = [fh_c*r_c/h_c; 0; h_c/r_c^2];
for n = 1:ObserverTest.Nagents,
    
    % Compute the Gaussian noise based on the prescribed covariances and
    % mean value (edit F.Oliva 07/10)
    GPSnoise = ObserverTest.GPSGaussianCovariance.*randn(6,1);
    if(ObserverTest.ZeroErrors == 1)
        GPSnoise = 0*GPSnoise;
    end
    
    if(ObserverTest.GPS_LVLH2ECI == 1)
        % inertial position and velocity into relative LVLH
        % frame back to ECI
        
        if(n == 1) %chief
            GPS_LVLH = GPSnoise;
        else
            GPS_LVLH = deputy_rel_LVLH(:,n-1) + GPSnoise; %GPS given wrt Chief's LHLV
        end
        myGPS = LVLH2ECI_V1_1(GPS_LVLH(1:3),coe(3), coe(5), coe(4) + coe(6)) + satellite_pos_ECI(:,1);
        myGPSpeed = LVLH2ECI_V1_1(GPS_LVLH(4:6) + cross(chief_omega_LVLH, GPS_LVLH(1:3)),coe(3), coe(5), coe(4) + coe(6)) + satellite_vel_ECI(:,1);
    else
        myGPS = GPSnoise(1:3)+ Agent(n).iner_ECI(ObserverTest.ExtractGPS,i+1);
        myGPSpeed = GPSnoise(4:6) + Agent(n).iner_ECI(4:6,i+1);
    end
    
    
    %IF GPS IS NOT READ CORRECTLY, IT IS REPLACED WITH THE
    %ESTIMATION plus noise...or past values....take care...???
    if(Agent(n).SuccessfullyReadGPS(end) == 0)
        %myGPS = Agent(n).xHatUKF(1:3,i+1) + 0.8*ObserverTest.GPSGaussianCovariance(1:3).*randn(3,1);
        %myGPSpeed = Agent(n).xHatUKF(4:6,i+1)+ 0.8*ObserverTest.GPSGaussianCovariance(4:6).*randn(3,1);
        myGPS = ObserverTest.ThetaLossGPSdata*(Agent(n).GPSopt(1:3,i) + Agent(n).GPSpeed(1:3,i)*ObserverTest.time_step) + ...
                        (1-ObserverTest.ThetaLossGPSdata)*Agent(n).xHatUKF(1:3,i+1);
        myGPSpeed =  ObserverTest.ThetaLossGPSdata*Agent(n).GPSpeed(1:3,i) + (1-ObserverTest.ThetaLossGPSdata)*Agent(n).xHatUKF(4:6,i+1);
    end
    
    
    %PROJECTION OF THE GPS (for now only position)
    %if(i>ObserverTest.UWBOptimizationNoBeforeThan && Agent(n).MismatchAprioriDistances(i) < ObserverTest.AprioriMismatch4Projection)
    if(i>ObserverTest.UWBOptimizationNoBeforeThan)% && Agent(n).SuccessfullyReadGPS(end)==1)
        dynamicSphere = ObserverTest.GPSprojectionOnAprioriSphereRadius*max(1,Agent(n).MismatchAprioriDistances(i) / ObserverTest.AprioriMismatch4Projection);
        if(norm(Agent(n).xHatUKF(1:3,i+1)-myGPS) > dynamicSphere)     %too far!
            franorm = norm(Agent(n).xHatUKF(1:3,i+1)-myGPS)/dynamicSphere;
            myGPS = Agent(n).xHatUKF(1:3,i+1)+(myGPS-Agent(n).xHatUKF(1:3,i+1))/franorm;
        end
    end
    

    
    Agent(n).GPS(:,i+1) = myGPS;
    Agent(n).GPSpeed(:,i+1) = myGPSpeed;
    for jj = 1:3,
        Agent(n).GPSfilter([1 2] +(jj-1)*2,i+1) =  ObserverTest.FilterAposition*Agent(n).GPSfilter( [1 2] +(jj-1)*2,i) +ObserverTest.FilterBposition*Agent(n).GPS(jj,i);
        Agent(n).GPSfiltered(jj,i+1) = ObserverTest.FilterDposition*myGPS(jj) + ObserverTest.FilterCposition*Agent(n).GPSfilter([1 2] +(jj-1)*2,i+1);
    end
    
    %calcolo della pseudoderivata sul segnale filtrato del
    %GPS
    IterativePseudoDerivativeCodeAgentGPS(i+1,n); %
    if(ObserverTest.AddDerivativeEstimation == 2) %GPS speed meseasured by GPS....!
        Agent(n).GPSderivative(:, i+1) = Agent(n).GPSpeed(:,i+1);
    end
    
    %specific constant bias of the agent
    if(ObserverTest.GPSoffset == 1 && ObserverTest.ZeroErrors == 0)
        myGPS = myGPS + Agent(n).GPSoffset;
    end
    
    Agent(n).GPSopt(:,i+1) = myGPS;
    Agent(n).GPS(:,i+1) = myGPS;
    Agent(n).GPSpeedOpt(:,i+1) = Agent(n).GPSpeed(:,i+1);
    
    %Data collection with a priori estimates of each agent
    %exchanged over the network
    for jj=1:ObserverTest.Nagents
        Agent(n).GotData(1:end-1,(jj-1)*3+1:(jj-1)*3+3) =  Agent(n).GotData(2:end,(jj-1)*3+1:(jj-1)*3+3);
        if(ObserverTest.TrasmittedTruePositions == 1)
            Agent(n).GotData(end,(jj-1)*3+1:(jj-1)*3+3) =  Agent(jj).iner_ECI(1:3,i+1);
        else
            Agent(n).GotData(end,(jj-1)*3+1:(jj-1)*3+3) =  Agent(jj).xHatUKF(1:3,i+1);
        end
    end
    Agent(n).MismatchAprioriDistances(i+1) = sum(abs(ObserverTest.MeasuredDistances(n,:)-ObserverTest.AprioriAdjacencyMatrix(n,:)));
    Agent(n).ErrorAprioriDistances(i+1) = sum(abs(ObserverTest.AdjacencyMatrix(n,:)-ObserverTest.AprioriAdjacencyMatrix(n,:)));
    
end
