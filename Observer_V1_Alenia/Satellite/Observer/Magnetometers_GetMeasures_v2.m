function [MagTemp,MagTemp2,MagECI] = Magnetometers_GetMeasures_v2(DynOpt,Agent_iner_eci,Agent_quaternion,Dtheta)

% This function provides the magnitude and gyros measurements if observerRequest==0 (data need to be the true one and not the estimated ones)
% adding estimated sensor noise, whereas if observerRequest == 1 it evaluates the shyntetic measures given the estimated vectors
% Dtheta: is the 3 components vector of roll, pitch, yaw in radiants that
% identifies the relative rotation between the first and the second
% magnetometer

if DynOpt.integration_pos == 1
    myutc = [2019 12 15 10 20 36]; %CHANGE THIS...??
    LatLongAlt = eci2lla(Agent_iner_eci*1E3,myutc); %converto from ECI to latitude, longitude,  altitude
    [mag_field_vector,~,~,~,~] = igrfmagm(max(1000,min(LatLongAlt(3),6E5)),LatLongAlt(1),LatLongAlt(2),decyear(2019,12,15),13); %mag_field_vector is in nanotesla, by IGRF11-12
    DynOpt.mag_field_vector = mag_field_vector;
else
    mag_field_vector = DynOpt.mag_field_vector;
end
q_ECI2Body =  Agent_quaternion; 
wrap_mag = (mag_field_vector'/(norm(mag_field_vector)));

R_ECI2Body = quat2dcm(q_ECI2Body);
% R_ECI2Body = eul2rotm(quat2eul(q_ECI2Body)-wrap_mag');

MagTemp = (R_ECI2Body)*(wrap_mag);

% 2nd Magnetometer, evaluate magnetic measurements 
R_ECI2Body_2 = eul2rotm(Dtheta)*R_ECI2Body ; %rotation matrix accounting the roll-pitch-yaw rotation that allows to pass from sensor 1 to sensor 2

%%% ADD NOISE BETWEEN SENSORS %%%
tempangle = randn(1,3)*DynOpt.EulerAngleNoiseOnMag;
R_ECI2Body_2 = eul2rotm(tempangle)*R_ECI2Body_2;

MagTemp2 = (R_ECI2Body_2)*(wrap_mag); %same bias for now....

%%% ECI magnetometer %%%
MagECI = wrap_mag;

end
                                                       

                                                       
                                                       
               