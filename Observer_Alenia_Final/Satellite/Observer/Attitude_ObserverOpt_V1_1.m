function  Attitude_ObserverOpt_V1_1(satellites_attitude, params,i,tspan_att,k)
%Attitude_Obsever_V2_1: estimates the i-th satellite attitude using a UKF
% joining the papers "A Simple Attitude Unscented Kalman Filter: Theory and
% Evaluation in a Magnetometer-Only Spacecraft Scenario" by MURTY S. CHALLA et al., Digital Object Identifier 10.1109/ACCESS.2016.2559445
% and the paper "Inexpensive CubeSat Attitude Estimation Using Quaternions
% and Unscented Kalman Filtering" by Vinther, Kasper; Fuglsang Jensen, Kasper; Larsen, Jesper Abildgaard; Wisniewski, Rafal
% Published in: Automatic Control in Aerospace
%
% INPUT:
% satellites_attitude: Array (7*N x 1), with N equal to the number of satellites. It is built in the following way: [quaternion1; omega1; quaternion2; omega2; etc...
%                    .... ], where the first component of the quaternion is    the scalar component
% deputy_rel_LVLH: deputy rLVLH coordinates with respect to chief
% params:  object defined by La Sapienza containing general costants,  simulation parameters, control torques, desired references...
% i: previous time step index (the estimates are computed at time t = (i+1)*time_step), IT IS ASSUMED THAT time_step_att = time_step  (position).
% tspan: Time interval used in the integration (prediction step) = [i*time_step, (i+1)*time_step]
% k: k-th satellites
%
%OUTPUT: the global variables Agent ObserverTest are updated with the
%attitude estimation

global Agent ObserverTest

%Evaluate the actual attitude estimation integrating one step ahead
%i+1
X =  rk4_V1_1_decentralized(@AttitudeDynamicsOpt_V1_1_decentralized, tspan_att, Agent(k).attitude_xHatUKF(:,i), params,k);
Agent(k).attitude_xHatUKF(:,i+1) =X(:,end);

%GET MEASURES, steps 2.1-2.7
%constant (or not) magneto and gyro bias
Agent(k).magnetoBias(:,i+1) = ObserverTest.MagnetoBias;
Agent(k).gyroBias(:,i+1) = ObserverTest.GyroBias;
%get synthetic measurements (true)
[MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(Agent(k).iner_ECI(1:3,i+1)',Agent(k).attitude(1:4,i+1)',...
    Agent(k).attitude(5:7,i+1)',Agent(k).magnetoBias(:,i+1)',Agent(k).gyroBias(:,i+1)',0); %????
Agent(k).magneto(:,i+1) = MagTemp;
Agent(k).gyros(:,i+1) = GyrosTemp;

ObserverTest.attitudeTimeIndex = i+1;
ObserverTest.Kth_agent = k;

if(ObserverTest.attitudeTimeIndex > (ObserverTest.attitudeBackSamplesWindow-1)*(ObserverTest.attitudeBackIntersample) )%enough data point to optimize bak in the data window
    ObserverTest.attitudeOptimizationCounter = 0;
    [NewState, Jattitude] = fminunc(@attitudeObsOptimization_V1_1,Agent(k).attitude_xHatUKF(:,i+1),ObserverTest.optionsObs);    
    disp(['Attitude, drone  ' num2str(k) ': Jstart = ' num2str(ObserverTest.attitudeOptimizationCyclesInitial) ...
            ', Jfinal/Jstart = ' num2str(ObserverTest.attitudeOptimizationCyclesEnd/ObserverTest.attitudeOptimizationCyclesInitial)]);
        %Agent(k).iner_ECI(ExtractGPS,i+1);
    %step 2.12
    Agent(k).attitude_xHatUKF(:,i+1) = NewState;
end


end

