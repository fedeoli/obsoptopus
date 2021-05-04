function  Attitude_ObserverOpt_V2_EKF(map, params,i,tspan_att,k)
%Attitude_Obsever_V2_1: estimates the i-th satellite attitude using MAHONY

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

% global vars
global Agent ObserverTest

%GET MEASURES, steps 2.1-2.7
%constant (or not) magneto and gyro bias
Agent(k).magnetoBias(:,i+1) = ObserverTest.MagnetoBias;
Agent(k).gyroBias(:,i+1) = ObserverTest.GyroBias;

%get synthetic measurements (true)
[MagTemp_body1, MagTemp_body2, GyrosTemp_body] = AttitudeObserver_GetMeasures_v2_2(Agent(k).iner_ECI(1:3,i+1)',Agent(k).attitude(1:4,i+1)',...
    Agent(k).attitude(5:7,i+1)',Agent(k).magnetoBias(:,i+1)',Agent(k).gyroBias(:,i+1)',ObserverTest.RPYbetweenMagSensors,0); 
Agent(k).magneto(:,i+1) = MagTemp_body1;
Agent(k).magneto2(:,i+1) = MagTemp_body2;
Agent(k).gyros(:,i+1) = GyrosTemp_body;

%%% Handling a priori attitude coordinates
% get quaternion evolution
q_true = Agent(k).attitude(1:4,i+1)'; 
% get angular velocity evolution
omega_true = Agent(k).attitude(5:7,i+1)';

%%% Apriori Estimated Quaternion Dynamics as in Mahony 10.1109/TAC.2008.923738
xhat =  rk4_V1_1_decentralized(@AttitudeDynamics_V1_2_decentralized, tspan_att, Agent(k).attitude_xHatUKF(:,i), params,k);

% get quaternion evolution
q_hat = xhat(1:4,end)'; 
% get angular velocity evolution
omega_hat = xhat(5:7,end)';

% quaternion normalization
q_ECI2Body = q_true/norm(q_true);
q_ECI2Body_hat = q_hat/norm(q_hat);

% get attitude values
attitude_ECI = q_ECI2Body(2:4);
attitude_ECI_est = q_ECI2Body_hat(2:4);

% reference transformation and normalization (ref 5.3)
R_ECI2Body = quat2dcm(q_ECI2Body); 
Agent(k).R_ECI2Body = R_ECI2Body;

% transformation on magnetometer measures
MagTemp_I1 = R_ECI2Body'*MagTemp_body1;
if ObserverTest.nMagneto == 2
    MagTemp_I2 = R_ECI2Body'*MagTemp_body2;
end

% state and oberved state
xhat = [q_ECI2Body_hat, omega_hat];
xhat_past = Agent(k).attitude_xHatUKF(1:7,i)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ObserverTest.Sun SENSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ObserverTest.Sun == 1
    [hss_I, hss_body] = Sun_measure(k,params);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numeric evaluaion of the output mapping
BI_1 = params.BI1_sym;
if ObserverTest.nMagneto == 2
    BI_2 = params.BI2_sym;
end
if ObserverTest.Sun == 1
    PSun = params.Psun;
end

% symlinks
f = map.f;
h = map.h;
A = map.A; % dynamic and transmission jacobian matrices 
H = map.H;
omega_Body2ECI_Body = map.omega_Body2ECI_Body;
q_sym = map.q_sym;

% output mapping - computation
h_num = subs(h, BI_1, MagTemp_I1');
if ObserverTest.nMagneto == 2
    h_num = subs(h_num, BI_2, MagTemp_I2');
end
if ObserverTest.Sun == 1
    h_num = subs(h_num, PSun, hss_I');
end
h_num = subs(h_num, omega_Body2ECI_Body, omega_hat);
h_num = subs(h_num, q_sym, q_ECI2Body_hat);
h_num = double(h_num);

% measurement assignment (ref 5.4)
if ObserverTest.Sun == 0
    if ObserverTest.nMagneto == 1 
        z = [MagTemp_body1' GyrosTemp_body'];
    elseif ObserverTest.nMagneto == 2
        z = [MagTemp_body1' MagTemp_body2' GyrosTemp_body'];
    end
elseif ObserverTest.Sun == 1
    if ObserverTest.nMagneto == 1 
        z = [MagTemp_body1' hss_body' GyrosTemp_body'];
    elseif ObserverTest.nMagneto == 2
        z = [MagTemp_body1' MagTemp_body2' hss_body' GyrosTemp_body'];
    else
        z = [hss_body' GyrosTemp_body'];
    end
end

% numeric evaluation (ref 5.5)
z_hat = double(h_num'); 

% Params storage procedure
params.R_ECI2Body = R_ECI2Body;
params.omega_Body2ECI_Body = omega_Body2ECI_Body;
params.q_ECI2Body = q_ECI2Body;
params.attitude_ECI = attitude_ECI;
params.attitude_ECI_est = attitude_ECI_est;
params.P = Agent(k).P_att;
params.Q = Agent(k).Q_att;
params.R = Agent(k).R_att;
params.X = map.X;
params.A = A;
params.H = H;
params.ObserverTest = ObserverTest;
params.tspan_att = tspan_att;

%%%%% EKF %%%%%
[x_hat_next, Pnext] = ekf_V5(xhat_past,xhat,z,z_hat,params,ObserverTest.nMagneto,ObserverTest.Sun);

% update cvariance matrix
Agent(k).P_att = Pnext;

% update estimation state + measures
xhat(1:4) = x_hat_next(1:4);
xhat(5:7) = x_hat_next(5:7);

% estimation assignment
Agent(k).attitude_xHatUKF(:,i+1) = [xhat'; zeros(3,1); zeros(3,1)];
end

