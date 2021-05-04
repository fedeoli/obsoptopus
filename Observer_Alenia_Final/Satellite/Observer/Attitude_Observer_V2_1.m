function  Attitude_Observer_V2_1(satellites_attitude, params,i,tspan_att,k)
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


%% ATTITUDE UKF
%[W0m,W0c,Wim,CHI] = AttitudeSigmaPoints_V1_1(Agent(k).attitude_xHatUKF(:,i),Agent(k).attitudeP);
%Agent(k).attitudeP = eye(length(Agent(k).attitudeP))*1E-6;%???TEST LEAVE IT!!!
[W0m,W0c,Wim,CHI] = AttitudeSigmaPoints_V1_1(Agent(k).attitude_xHatUKF(:,i),Agent(k).attitudeP);

%Sigma Points integration  (a-priori sigma)
Lukf = length(Agent(k).attitude_xHatUKF(:,i))-1; %the first component of the quaternion is constrained by the others three that's way -1
%tspan_att = time(i): time_step_att: time(i) + time_step;
CHI_meno = zeros(Lukf+1,2*Lukf+1); %number of rows is Lukf + 1 because of the first element of the quaternion!
for j=1:2*Lukf+1,
    %CHI(:,j) = [Agent(k).attitude(:,i+1);zeros(6,1)];%CHI(:,1);%???TEST LEAVE IT!!!
    X_Attitude = rk4_V1_1_decentralized(@AttitudeDynamics_V1_2_decentralized, tspan_att, CHI(:,j), params,k);
    CHI_meno(:,j) = X_Attitude(1:end,end);
end
%A priori estimate, step 1.5 of paper
WIM = [W0m,Wim*ones(1,2*Lukf)]';
attitude_xHatMeno = CHI_meno*WIM; %x_k|k-1 %but for quaternion things should be done differently as stated in
%attitude_xHatMeno(1:4,1) = quatnormalize(attitude_xHatMeno(1:4,1)');

%the paper A Simple Attitude Unscented Kalman Filter: Theory and Evaluation in a Magnetometer-Only Spacecraft Scenario
%MURTY S. CHALLA et al., Digital Object Identifier 10.1109/ACCESS.2016.2559445
chi_meno_invquaternion = quatinv(CHI_meno(1:4,1)');
average_angular_deviation_half =  zeros(1,3); %eq. (57) MURTY S. CHALLA
for j=2:2*Lukf+1,
    temp = quatmultiply(chi_meno_invquaternion,CHI_meno(1:4,j)'); %eq. (56) MURTY S. CHALLA
    average_angular_deviation_half = average_angular_deviation_half + WIM(j)*temp(2:4);
end
attitude_xHatMeno(1:4,1) = quatnormalize(quatmultiply(CHI_meno(1:4,1)',[1, average_angular_deviation_half])); %eq. (58) MURTY S. CHALLA
q_apriori = attitude_xHatMeno(1:4,1);
omega_apriori = attitude_xHatMeno(5:7,1);
MagBias_apriori = attitude_xHatMeno(8:10,1);
GyrosBias_apriori = attitude_xHatMeno(11:13,1);

%Vinther
%                 deltaCHImeno = 0*CHI_meno; %we will not use the first component that is related to the first element of the quaternion
%                 for j=1:2*Lukf+1,%step 1.6 but
%                     deltaCHImeno(:,j) = [quatmultiply(CHI_meno(1:4,j)', quatinv(attitude_xHatMeno(1:4,1)') )...
%                                                               (CHI_meno(5:end,j) - attitude_xHatMeno(5:end,1))']';
%                 end
%                 deltaCHImeno = deltaCHImeno(2:end,:); %we clear off the first scalar component of the quaternion
%                 delta_x_meno = deltaCHImeno*WIM; %step 1.7

%as in CHALLA, eq. (59)
xi = zeros(Lukf,2*Lukf+1);
for j=1:2*Lukf+1,
    temp =  quatmultiply(quatinv(q_apriori'),CHI_meno(1:4,j)');
    xi(1:3,j) = 2*temp(2:4)'; %2*temp(2:4)'; ????
    xi(4:end,j) = CHI_meno(5:end,j)-attitude_xHatMeno(5:end);
end

%a priori covariance matrix
%Vinther
%                 attitude_Pmeno = ObserverTest.AttitudeQ;
%                 for j=1:2*Lukf+1,%step 1.8, we need to skip the first component
%                     attitude_Pmeno = attitude_Pmeno + WIM(j)*(deltaCHImeno(:,j)-delta_x_meno)*(deltaCHImeno(:,j)-delta_x_meno)';
%                 end
%as in CHALLA, eq. (60) but with Q since a different estended state as been considered, as in the Vinther paper
attitude_Pmeno = zeros(Lukf);
for j=1:2*Lukf+1,
    attitude_Pmeno = attitude_Pmeno + WIM(j)*xi(:,j)*xi(:,j)';
end

%GET MEASURES, steps 2.1-2.7
%constant (or not) magneto and gyro bias
Agent(k).magnetoBias(:,i+1) = ObserverTest.MagnetoBias;
Agent(k).gyroBias(:,i+1) = ObserverTest.GyroBias;
 
%get synthetic measurements (true)
[MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(Agent(k).iner_ECI(1:3,i+1)',Agent(k).attitude(1:4,i+1)',...
    Agent(k).attitude(5:7,i+1)',Agent(k).magnetoBias(:,i+1)',Agent(k).gyroBias(:,i+1)',0); %????
Agent(k).magneto(:,i+1) = MagTemp;
Agent(k).gyros(:,i+1) = GyrosTemp;
z = [MagTemp;GyrosTemp]; %vector of measurements
%z = [Agent(k).attitude(1:7,i+1)];%???? direct measure

%Perform synthetic measures from CHI points generating the
%corresponding output, %step 2.8
Z = zeros(6,2*Lukf+1);
%Z = zeros(6+1,2*Lukf+1);%???? direct measure
for j=1:2*Lukf+1, %eq. (27)
    %measurements associated to the sigma points
    [ZMag,ZGyros] = AttitudeObserver_GetMeasures_v2_1(Agent(k).xHatUKF(1:3,i+1)',CHI_meno(1:4,j)',...
        CHI_meno(5:7,j)',0*CHI_meno(8:10,j)',0*CHI_meno(11:13,j)',1);%?????
    Z(:,j) = [ZMag;ZGyros];
    %Z(:,j) = CHI_meno(1:7,j);%???? direct measure
end
zHatMeno = Z*WIM; %eq. (28)

Pzz = ObserverTest.AttitudeR;%eq. (31)
%Pzz = blkdiag(ObserverTest.AttitudeR(1,1),ObserverTest.AttitudeR);%eq. (31), %???? direct measure
for j=1:2*Lukf+1,
    Pzz = Pzz  + WIM(j)*(Z(:,j)-zHatMeno)*(Z(:,j)-zHatMeno)';
end

Pxz = zeros(Lukf,6);
%Pxz = zeros(Lukf,6+1);%???? direct measure
for j=1:2*Lukf+1, %eq. (32)
    %Vinther
    %Pxz = Pxz  + WIM(j)*(CHI_meno(2:end,j)-attitude_xHatMeno(2:end))*(Z(:,j)-zHatMeno)';
    %as in CHALLA, eq. (60)
    Pxz = Pxz  + WIM(j)*xi(:,j)*(Z(:,j)-zHatMeno)';
end

%step 2.9, the Kalman gain
attitudeK = Pxz*pinv(Pzz);
%step 2.10
delta_attitude_xHatMeno = attitudeK*(z - zHatMeno);
delta_q_aposteriori = delta_attitude_xHatMeno(1:3);
delta_omega_aposteriori = delta_attitude_xHatMeno(4:6);
delta_MagBias_aposteriori = delta_attitude_xHatMeno(7:9);
delta_GyrosBias_aposteriori = delta_attitude_xHatMeno(10:12);
%step 2.11, expand quaternion
%q_aposteriori = quatmultiply([sqrt(1-delta_q_aposteriori'*delta_q_aposteriori) delta_q_aposteriori' ] , q_apriori');
%modified as in CHALLA (66)-(67)
delta_q_update = [1, delta_attitude_xHatMeno(1:3)'/2.0];
q_aposteriori =quatnormalize(quatmultiply(q_apriori', delta_q_update));

%step 2.12
Agent(k).attitude_xHatUKF(:,i+1) = [q_aposteriori'; omega_apriori+delta_omega_aposteriori;...
    MagBias_apriori+delta_MagBias_aposteriori;GyrosBias_apriori+delta_GyrosBias_aposteriori;];
%Agent(k).attitude_xHatUKF(:,i+1) = [q_aposteriori'; omega_apriori+delta_omega_aposteriori;...
%                                                     zeros(3,1);zeros(3,1)];%????LEVA!!!!!
%step 2.13
temp = Agent(k).attitudeP;
Agent(k).attitudeP = attitude_Pmeno - attitudeK*Pzz*attitudeK';
if(min(eig(Agent(k).attitudeP))<0)
    Agent(k).attitudeP = temp;
end


end

