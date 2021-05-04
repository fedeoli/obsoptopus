function [MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(Agent_iner_eci,Agent_quaternion,Agent_omega,...
                                                                                                                    Agent_MagBias, Agent_GyroBias,UKF)

% This function provides the magnitude and gyros measurements if UKF==0 (data need to be the true one and not the estimates)
% adding estimated sensor noise, 
% whereas if UKF == 1 it evaluates the shyntetic measures given the estimated vectors

% Agent(k_agent).iner_ECI(1:3) =  Agent_iner_eci
% Agent(k_agent).attitude(1:4,time..) = Agent_quaternion
% Agent(k_agent).attitude(5:7,time..) = Agent_omega

global ObserverTest 

MagBias = Agent_MagBias';
GyroBias = Agent_GyroBias';
EulerAngleNoiseOnMag = ObserverTest.EulerAngleNoiseOnMag;
ErrorMu = 1;
if(ObserverTest.AttitudeZeroErrors==1 && UKF == 0)
    ErrorMu = 0;
    MagBias = 0*MagBias;
    GyroBias = 0*GyroBias;
    EulerAngleNoiseOnMag = 0;
else
    if(UKF == 1)
        %MagBias = Agent_MagBias';%LEAVE THE ZERO!!! ????
        %GyroBias = Agent_GyroBias';%LEAVE THE ZERO!!! ????
        EulerAngleNoiseOnMag = 0;
    end
end

myutc = [2019 12 15 10 20 36]; %CHANGE THIS...??
LatLongAlt = eci2lla(Agent_iner_eci*1E3,myutc); %converto from ECI to latitude, longitude,  altitude
[mag_field_vector,hor_intensity,declinatioon,inclination,total_intensity] ...%IT IS SATURATED FOR MOST TRAJECTORY(ALTITUDE)!!!!
                 = igrfmagm(max(1000,min(LatLongAlt(3),6E5)),LatLongAlt(1),LatLongAlt(2),decyear(2019,12,15),12); %
                %mag_field_vector is in nanotesla, by IGRF11-12
mag_field_ECI = mag_field_vector;
q_ECI2Body =  Agent_quaternion; 

%ADDING NOISE directly on Euler angles

% if(UKF == 0) %noise is added only when a synthetic measurement is requested
%     [r,p,y] = quat2angle(q_ECI2Body);
%     temp = [r,p,y]'+randn(3,1)*EulerAngleNoiseOnMag;
%     q_ECI2Body = quatnormalize(angle2quat(temp(1),temp(2),temp(3)));
% end

R_ECI2Body = quat2dcm(q_ECI2Body) ;

%Magnetic readigs in the Body frame as would be done by a
%real sensor mounted on the cubesat, so in magneto there is
%the measure performed by the sensor (sensor simulation is done with "correct" data not the estimated ones).

MagTemp = R_ECI2Body*mag_field_vector'/(norm(mag_field_vector)) + MagBias;
%add noise with conversion [azimuth,elevation,r] = cart2sph(X,Y,Z) and then back [x,y,z] =sph2cart(azimuth,elevation,r).
if(UKF == 0) 
    %[azimuth,elevation,r] = cart2sph(MagTemp(1),MagTemp(2),MagTemp(3));
    %[x,y,z] = sph2cart(azimuth  +ErrorMu*(ObserverTest.headingMagDisplacement + randn(1)*ObserverTest.MagSigmaAzimut),...
     %                                                      elevation+ErrorMu*(ObserverTest.tiltMagDisplacement + randn(1)*ObserverTest.MagSigmaElevation),r);
    %MagTemp =[x,y,z]';
end

%Angular velocities measures are simpler
%Agent(k_agent).attitude(5:7,time_i); is the ^c omega(t) of the refererred paper
GyrosTemp = Agent_omega'  + GyroBias + ErrorMu*ObserverTest.GyroGaussianCovariance'.*randn(3,1);
if(UKF == 0)
   %[azimuth,elevation,r] = cart2sph(GyrosTemp(1),GyrosTemp(2),GyrosTemp(3));
   %[x,y,z] = sph2cart(azimuth  +ErrorMu*ObserverTest.headingGyroDisplacement,...
   %                                                        elevation+ErrorMu*ObserverTest.tiltGyroDisplacement ,r);
   %GyrosTemp =[x,y,z]';
end
                                                       

                                                       
                                                       
               