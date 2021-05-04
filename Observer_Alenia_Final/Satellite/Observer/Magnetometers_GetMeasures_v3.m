function [MagTemp,MagTemp2,DynOpt] = Magnetometers_GetMeasures_v3(DynOpt,pos,Agent_quaternion,Dtheta)

% This function provides the magnitude and gyros measurements if observerRequest==0 (data need to be the true one and not the estimated ones)
% adding estimated sensor noise, whereas if observerRequest == 1 it evaluates the shyntetic measures given the estimated vectors
% Dtheta: is the 3 components vector of roll, pitch, yaw in radiants that
% identifies the relative rotation between the first and the second
% magnetometer

%%% convert quaternion %%%
% Agent_quaternion = ConvertQuat_V2_1(Agent_quaternion, 'ScalarTo4');

%%% get magnetic field %%%
DynOpt.mag_field_vector = DynOpt.mag_field_story(:,pos);
mag_field_vector = DynOpt.mag_field_story(:,pos);

% quaternion and magnetic field
q_ECI2Body =  Agent_quaternion; 
wrap_mag = (mag_field_vector/(norm(mag_field_vector)));

% first magnetometer
% R_ECI2Body = RotationConversion_V2_1('QtoDCM', q_ECI2Body);
R_ECI2Body = quat2dcm(q_ECI2Body);

MagTemp = (R_ECI2Body)*(wrap_mag);

% 2nd Magnetometer, evaluate magnetic measurements 
% Rtheta = RotationConversion_V2_1('EA321toDCM', Dtheta');
Rtheta = eul2rotm(transpose(Dtheta));
R_ECI2Body_2 = Rtheta*R_ECI2Body ; 
MagTemp2 = (R_ECI2Body_2)*(wrap_mag); 

%%% ECI magnetometer measured%%%
DynOpt.MagECI = (wrap_mag);
DynOpt.Rmag_1 = R_ECI2Body;
DynOpt.Rmag_2 = R_ECI2Body_2;

end
                                                       

                                                       
                                                       
               