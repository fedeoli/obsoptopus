%% EKF PROJECT 
% Author: Daniele Giorgi, Federico Oliva
% Date: 6/09/2020
%
% POSSIBILI UPDATE/PROBLEMI
%
% 4) Il modello Albedo sembra funzionare bene, per testarlo però con la
% configurazione attuale del satellite_iner_ECI parte dall'eclissi, quindi
% mentre se si usa come coordinate sferiche del sole quelle date dalla tesi
% di riferimento:
%
%sunsph = [0 1.17 CONST.ua]; 
%
%allora si può apprezzare l'effetto Albedo che sporca effettivamente le
%misure. Altrimenti bisogna simulare perlomeno un'orbita completa.
%Bisogna inoltre aggiungere i plot per differenti altitudini e o longitudini,
% ma può essere fatto facilmente attraverso la funzione ecef2lla(sat_ecef),
%in quanto è stato tutto convertito in Earth-Centered-Earth_Fixed coordinates. 
%% Init Section


% Clear workspace
clc
clear all
close all

test = 2


    if test == 0
        nMagneto = 1;
        Sun = 1;
    else
        if test == 1
            nMagneto = 2;
            Sun = 0;
        else
            nMagneto = 2;
            Sun = 1;
        end
    end
    

% Symbolic definitions
fprintf('Setting Observer parameters\n');
syms Ixx Iyy Izz;       % System Inertia
syms wx wy wz;          % System angular velocity
syms q1 q2 q3 q0;       % Quaternion
syms Bx1 By1 Bz1 Bx2 By2 Bz2;          % Magnetometer measures

%Set params ObserverTest
addpath(genpath([pwd 'Observer']));
addpath(genpath([pwd 'Lib']));
addpath(genpath([pwd 'AlbedoToolbox-1.0']));
addpath(genpath(pwd));
Input1;
OutputInitialization_V2_2;
params.SatellitesCoordinates = satellites_iner_ECI;
params.u = u;
Observer_Setup_Parameter_Save_V1_1;



%%%%% ORBITAL PROPAGATOR  %%%%%%
% Initial satellite position 
Chi = satellites_iner_ECI;
satellites_iner_ECI = satellites_iner_ECI(1:3,1)';
satellites_iner_ECI_est = satellites_iner_ECI;

% flags
% nMagneto = 2; % number of Magnetometers (max. 2)
% Sun = 0; % 0: no Sun Sensor; 1: with Sun Sensor
plot_section = 0;
randerr = 0;

% error settings
max_deg_var = -20;   % deg
max_vel_var = 10;   % deg/s

% Simulation data
nIter = 120;
tspan = [0 1]; 

% Global vars definition
global ObserverTest

%% INITIAL CONDITIONS
%  x_m: previous estimated state [quaternion, omegas, bias_mag, bias_gyro],  the first component of the quaternion is the scalar part
%  P_x: previous covariance matrices

% Kinematic and measures init 

% angular velocity measured by gyro - Body frame
omega_Body2ECI_Body = [wx wy wz];
w_body = omega_Body2ECI_Body;

% quaternion in Body frame
q_ECI2Body = [q0 q1 q2 q3];
q_sym = [q0 q1 q2 q3];
% attitude from quaternion - vector elements
attitude_ECI = q_ECI2Body(2:4);
% junk assignment (?)
qin = q_ECI2Body;

% Magneto measures - Inertial frame
Magneto1 = [Bx1 By1 Bz1];
BI_1 = Magneto1;
params.BI1_sym = BI_1;

if nMagneto == 2
    Magneto2 = [Bx2 By2 Bz2];
    BI_2 = Magneto2;
    params.BI2_sym = BI_2;
end

% Dynamics init - inertia and mass
In = [Ixx, Iyy, Izz];

% Cosines matrix from quaternion
dcm = sym('dcm', [3,3]);

% MATLAB VERSION
dcm(1,1,:) = qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2;
dcm(1,2,:) = 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4));  
dcm(1,3,:) = 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3));
dcm(2,1,:) = 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4));
dcm(2,2,:) = qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2;
dcm(2,3,:) = 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2));
dcm(3,1,:) = 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3));
dcm(3,2,:) = 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2));
dcm(3,3,:) = qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2;

% Magneto measures - Body frame
B_body_1 = dcm*transpose(BI_1);

if nMagneto == 2
    B_body_2 = dcm*transpose(BI_2);
end

% Angular velocity - Inertial frame
wI = transpose(dcm)*transpose(omega_Body2ECI_Body);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Sun == 1
    Sun_init
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vett. di stato
X = [q_ECI2Body, omega_Body2ECI_Body]; 

if Sun == 0
    if nMagneto == 1 
        h = [B_body_1; transpose(w_body)];
    elseif nMagneto == 2
        h = [B_body_1 ; B_body_2; transpose(w_body)];
    end
elseif Sun == 1
    if nMagneto == 1 
        h = [B_body_1; hss; transpose(w_body)];
    elseif nMagneto == 2
        h = [B_body_1 ; B_body_2; hss; transpose(w_body)];
    else
        h = [hss; transpose(w_body)];
    end
end

%number of state
n=size(X,2);     
%number f observed values [attitude, gyros];
m =size(h,1);    

%standard deviation (std) of process f
q = 1*1E-6;    
%std of measurement h
r = 3*1E-3;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%q = 0;
%r = 0;
%b_mag = 0;
%b_gyro = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAFE VERSION
% covariance matrix of process
Q=q*eye(n); 
% covariance matrix of measurement
R=r*eye(m); 
P = blkdiag(eye(4), eye(3))*1E-3;

% Initial condition 
xk = satellites_attitude(1:7,1)';

% Init attitude
q0 = quatnormalize(xk(1:4));
% q0_euler_rad = ConvertAttitude(q0', 'quaternion', '321')';
q0_euler_rad = quat2eul(q0, 'ZYX')';
w0 = xk(5:7);

% Noise definition - deg
amp_noise_attitude = max_deg_var*pi/180;        % rad
amp_noise_velocity = max_vel_var*pi/180;        % rad/s

% type of error
if randerr == 1
    qhat_euler_rad = q0_euler_rad + amp_noise_attitude*randn(3,1);
    qhat = quatnormalize(eul2quat(qhat_euler_rad','ZYX'));
    w_hat = w0 + amp_noise_velocity*randn(1,3);
else
    qhat_euler_rad = q0_euler_rad + amp_noise_attitude;
    qhat = quatnormalize(eul2quat(qhat_euler_rad','ZYX'));
    w_hat = w0 + amp_noise_velocity;
end
% error init

xhat = [qhat w_hat];

% error definition
err_euler_deg = (q0_euler_rad - qhat_euler_rad)*180/pi;
err_w = w0 - w_hat;


% measurement initialization
% Sensor Measurements, true(synthetic) and estimate
% displacement between the 2 magnetometers - euler angle
% riga 130 di Observer_Setup_Parameter_V1
Dtheta = ObserverTest.RPYbetweenMagSensors;
% general measure
[MagTemp_body1,MagTemp_body2,GyrosTemp_body] = ...
    AttitudeObserver_GetMeasures_v2_2(satellites_iner_ECI, q0, ...
    w0,zeros(3,1)',zeros(3,1)',Dtheta, 1);

% reference transformation and normalization
R_ECI2Body = quat2dcm(q0);
MagTemp_I1 = R_ECI2Body'*MagTemp_body1;

if nMagneto == 2
    MagTemp_I2 = R_ECI2Body'*MagTemp_body2;
end

% measurement assignment
if Sun == 1
    hss_body = I0*[1 1 1]';
    hss_I = R_ECI2Body'*hss_body;
end

if Sun == 0
    if nMagneto == 1 
        z = [MagTemp_body1' GyrosTemp_body'];
    elseif nMagneto == 2
        z = [MagTemp_body1' MagTemp_body2' GyrosTemp_body'];
    end
elseif Sun == 1
    if nMagneto == 1 
        z = [MagTemp_body1' hss_body' GyrosTemp_body'];
    elseif nMagneto == 2
        z = [MagTemp_body1' MagTemp_body2' hss_body' GyrosTemp_body'];
    else
        z = [hss_body' GyrosTemp_body'];
    end
end
z_hat = z + r*randn(1,length(z));

%% Dynamics setup section
% satellite equation
f = attitude_Kin_eqs(omega_Body2ECI_Body, q_ECI2Body.', In, params, dcm); 

% Numerical substitution
% f = subs(f, M, [0.89, .87, .61]);
f = subs(f, In, [params.sat(1).I(1,1) params.sat(1).I(2,2) params.sat(1).I(3,3)]);

% linearization of satellite equations
A = jacobian(f,X);    
H = jacobian(h,X);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = [A1*transpose(Magneto) A2*transpose(Magneto) A3*transpose(Magneto) ...
%     A4*transpose(Magneto) zeros(3,3); zeros(3,7)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage array init

x_hatV = zeros(nIter+1,n);      
zV = zeros(nIter+1,m);
z_hatV = zeros(nIter+1,m);
xV = zeros(nIter+1,n);
ECI_store = zeros(nIter+1,3);
MAG_store = zeros(nIter+1,3);
Idiff_store = zeros(nIter+1,3);

% store EKF results
x_hatV(1,:) = xhat;
xV(1,:) = xk;
zV(1,:) = z;
z_hatV(1,:) = z_hat;

% store ECI and magnetometer and Sun Sensor
ECI_store(1,:) = satellites_iner_ECI;
MAG_store(1,:) = MagTemp_I1;
Idiff_store(1,:) = [0;0;0];

 % Error array init
 
        q_meas = zeros(nIter+1,4);
        q_est = zeros(nIter+1,4);
        q_meas_Euler = zeros(nIter+1,3);
        q_est_Euler = zeros(nIter+1,3);
        sun_error = zeros(nIter+1, 3);
        delq = zeros((nIter+1),4);
        delq_Euler = zeros((nIter+1),3);      


% CALCOLO ERRORE DI QUATERNIONE - init
q_meas(1,:) = xV(1,1:4);
q_meas(1,:) = transpose(quatnormalize(q_meas(1,:)));
% q_meas_Euler(1,:) = ConvertAttitude(q_meas(1,:)', 'quaternion', '321')';
q_meas_Euler(1,:) = quat2eul(q_meas(1,:), 'ZYX')';

q_est(1,:) = x_hatV(1,1:4);
q_est(1,:) = transpose(quatnormalize(q_est(1,:)));
% q_est_Euler(1,:) = ConvertAttitude(q_est(1,:)', 'quaternion', '321')';
q_est_Euler(1,:) = quat2eul(q_est(1,:), 'ZYX')';

delq_Euler(1,:) = q_meas_Euler(1,:)-q_est_Euler(1,:);
delq(1,:) = quatmultiply(quatconj(q_meas(1,:)), q_est(1,:));


%% DYNAMICS INTEGRATION SECTION
fprintf('Dynamics integration\n');
tic
for k=1:nIter
    
    % Info
    if mod(k,10) == 0
        fprintf('Iteration N: %d\n',k);
    end
    
    % Integrazione posizione del satellite
    Chi = rk4_V1_1(@InertialDynamicsIntegrator_V2_1, tspan, Chi(:,end), params);
    satellites_iner_ECI = Chi(1:3,end)';
    
    % Store past values
    x_past = xk;
    xhat_past = xhat;
    
    %%%%%%%%%%%%%%%%%%%% ATTITUDE PROPAGATION %%%%%%%%%%%%%%%%%%%%
    % Simulo evoluzione VERA dello stato
    % No need to change tspan as the system is shift time invariant
    x = rk4_V1_1(@AttitudeDynamics_V1_fede, tspan, xk, params)';

    % get quaternion evolution
    q_true = x(end,1:4); 
    % get angular velocity evolution
    omega_true = x(end,5:7);
    
    % Simulation of the system with added disturbances
    % Simulo evoluzione STIMATA dello stato
    xhat = rk4_V1_1(@AttitudeDynamics_V1_fede, tspan, xhat, params)';

    % get quaternion evolution
    q_hat = xhat(end,1:4); 
    % get angular velocity evolution
    omega_hat = xhat(end,5:7);
    
    % quaternion normalization
    q_ECI2Body = q_true/norm(q_true);
    q_ECI2Body_hat = q_hat/norm(q_hat);
    
    % get attitude values
    attitude_ECI = q_ECI2Body(2:4);
    attitude_ECI_est = q_ECI2Body_hat(2:4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sensor Measurements, true(synthetic) and estimate  
    % displacement between the 2 magnetometers - euler angle
    % riga 130 di Observer_Setup_Parameter_V1
    Dtheta = ObserverTest.RPYbetweenMagSensors;
    % general measure
    [MagTemp_body1,MagTemp_body2,GyrosTemp_body] = ...
        AttitudeObserver_GetMeasures_v2_2(satellites_iner_ECI, q_ECI2Body, ...
        omega_true,zeros(3,1)',zeros(3,1)',Dtheta, 0);
    
    % reference transformation and normalization
    R_ECI2Body = quat2dcm(q_ECI2Body); 
    MagTemp_I1 = R_ECI2Body'*MagTemp_body1;
    if nMagneto == 2
        MagTemp_I2 = R_ECI2Body'*MagTemp_body2;
    end
    
    % state and oberved state
    xk_model = [q_ECI2Body, omega_true];
    xk = [q_ECI2Body, GyrosTemp_body'];
    xhat = [q_ECI2Body_hat, omega_hat];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUN SENSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Sun == 1
        Sun_measure
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % numeric evaluation of system linearization
    H_num = subs(H, BI_1, MagTemp_I1');
    if nMagneto == 2
        H_num = subs(H_num, BI_2, MagTemp_I2');
    end
    if Sun == 1
        H_num = subs(H_num, Psun, hss_I');
    end
    H_num = subs(H_num, omega_Body2ECI_Body, omega_hat);
    H_num = subs(H_num, q_sym, q_ECI2Body_hat);
    H_num = double(H_num);
    
    % numeric evaluaion of the output mapping
    h_num = subs(h, BI_1, MagTemp_I1');
    if nMagneto == 2
        h_num = subs(h_num, BI_2, MagTemp_I2');
    end
    if Sun == 1
        h_num = subs(h_num, Psun, hss_I');
    end
    h_num = subs(h_num, omega_Body2ECI_Body, omega_hat);
    h_num = subs(h_num, q_sym, q_ECI2Body_hat);
    h_num = double(h_num);
    
    % measurement assignment
    if Sun == 0
        if nMagneto == 1 
            z = [MagTemp_body1' GyrosTemp_body'];
        elseif nMagneto == 2
            z = [MagTemp_body1' MagTemp_body2' GyrosTemp_body'];
        end
    elseif Sun == 1
        if nMagneto == 1 
            z = [MagTemp_body1' hss_body' GyrosTemp_body'];
        elseif nMagneto == 2
            z = [MagTemp_body1' MagTemp_body2' hss_body' GyrosTemp_body'];
        else
            z = [hss_body' GyrosTemp_body'];
        end
    end

    % estimation
%     z_hat = double(H_num*xhat')';
    z_hat = double(h_num');

    % Params storage procedure
    params.R_ECI2Body = R_ECI2Body;
    params.omega_Body2ECI_Body = omega_Body2ECI_Body;
    params.q_ECI2Body = q_ECI2Body;
    params.attitude_ECI = attitude_ECI;
    params.attitude_ECI_est = attitude_ECI_est;
    params.P = P;
    params.Q = Q;
    params.R = R;
    params.X = X;
    params.A = A;
    params.H = H;
    params.ObserverTest = ObserverTest;
    params.tspan = tspan;
    
    %%%%% EKF %%%%%
    [x_hat_next, Pnext] = ekf_V5_fede(xhat_past,xhat,z,z_hat,params,nMagneto,Sun);
    
    
    % store EKF results
    x_hatV(k+1,:) = x_hat_next;
    xV(k+1,:) = xk;
    zV(k+1,:) = z;
    z_hatV(k+1,:) = z_hat;
    ECI_store(k+1,:) = satellites_iner_ECI;
    MAG_store(k+1,:) = MagTemp_I1;
    if Sun == 1
        Idiff_store(k+1,:) = Idiff;
    end
    
    % update cvariance matrix
    P = Pnext;
    
    % update estimation state + measures
    xhat(1:4) = x_hat_next(1:4);
    xhat(5:7) = x_hat_next(5:7);
    
    % CALCOLO ERRORE DI QUATERNIONE
    q_meas(k+1,:) = xV(k+1,1:4);
    q_meas(k+1,:) = transpose(quatnormalize(q_meas(k+1,:)));
%     q_meas_Euler(k+1,:) = ConvertAttitude(q_meas(k+1,:)', 'quaternion', '321')';
    q_meas_Euler(k+1,:) = quat2eul(q_meas(k+1,:), 'ZYX')';
    
    q_est(k+1,:) = x_hatV(k+1,1:4);
    q_est(k+1,:) = transpose(quatnormalize(q_est(k+1,:)));
%     q_est_Euler(k+1,:) = ConvertAttitude(q_est(k+1,:)', 'quaternion', '321')';
    q_est_Euler(k+1,:) = quat2eul(q_est(k+1,:), 'ZYX')';
    
    delq_Euler(k+1,:) = q_meas_Euler(k+1,:)-q_est_Euler(k+1,:);
    delq(k+1,:) = quatmultiply(quatconj(q_meas(k+1,:)), q_est(k+1,:));
    
   

    %%%%% Angular Sun error %%%%% 
    if Sun == 1
        sun_error(k,:) = sqrt(mean((S_true*180/pi-S_hat*180/pi).^2));
    end
end
    name = sprintf('delq_Euler%d',test);
    save(name, 'delq_Euler');

   if test == 2
        plot_section = 1;
    end
    toc
    
    %% Plot section
    if plot_section == 1
        %                 figure
        %                 for k=1:3
        %                     %set(gca, 'YAxisLocation', 'origin') % plot results
        %                     subplot(3,1,k)
        %                     plot3(1:(nIter+1), q_meas_Euler(:,k)*180/pi, '-', 1:(nIter+1), q_est_Euler(:,k)*180/pi, '--');
        %                     grid on;
        %                     lab = ["Roll","Pitch","Yaw"];
        %                     title('Roll-Pitch-Yaw attitude estimation');
        %                     legend(lab(k),strcat(lab(k), ' estimated'));
        %                 end
        %
        lab = ["Roll","Pitch","Yaw"];

        delq_Euler0 = load('delq_Euler0.mat');
        delq_Euler1 = load('delq_Euler1.mat');
        delq_Euler2 = load('delq_Euler2.mat');

        figure
        
        for k=1:3
            subplot(3,1,k)
            plot(1:(nIter+1), abs(delq_Euler0.delq_Euler(1:nIter+1,k))*180/pi, 'r-', ...
                 1:(nIter+1), abs(delq_Euler1.delq_Euler(1:nIter+1,k))*180/pi, 'b-', ...
                 1:(nIter+1), abs(delq_Euler2.delq_Euler(1:nIter+1,k))*180/pi, 'g-');
            grid on;
            %quat = sprintf('q%d',k);
            %strcat(quat,' error'),
            title(lab(k)+ ' error');
            %legend(lab(k));
            legend('Magneto = 1, Sun = 1', ...
                'Magneto = 2, Sun = 0', ...
                'Magneto = 2, Sun = 1');
            %axis([20 120 -inf +inf])
            xlabel('time (s)');
            ylabel('error(deg)');

        end
        
        figure
        for k=1:3
            subplot(3,1,k)
            plot(1:(nIter+1), abs(delq(1:nIter+1,k)), 'r-', ...
                 1:(nIter+1), abs(delq((nIter+1):(2*nIter+1),k)), 'b-', ...
                 1:(nIter+1), abs(delq((2*nIter+1):(3*nIter+1),k)), 'g-');
            grid on;
            quat = sprintf('q%d',k+1);
            strcat(quat,' error');
            title('Attitude error');
            legend('Magneto = 1, Sun = 1', ...
                'Magneto = 2, Sun = 0', ...
                'Magneto = 2, Sun = 1');
            %legend(quat);
            %ylim([-2 8])
        end
% real trajectory
figure
grid on
hold on
title('Fleet trajectories')
xlabel('Km')
ylabel('Km')
zlabel('Km')
RF_flag = 0;
style_traj = 'o';
linewidth = 1;
% satellite
xreal = ECI_store(:,1);
yreal = ECI_store(:,2);
zreal = ECI_store(:,3);
plot3(xreal,yreal,zreal,style_traj,'MarkerFaceColor','b');

% magnetic field measurements - Magnetometer 1
figure
hold on
plot(1:nIter+1,MAG_store(:,1),'r');
plot(1:nIter+1,MAG_store(:,2),'g');
plot(1:nIter+1,MAG_store(:,3),'b');
legend('Mag_I_x','Mag_I_y','Mag_I_z');

% difference between model with albedo and without
figure
hold on
for i=2:nIter+1
plot(i,norm(Idiff_store(i,:)),'bo')
end

figure
for k=1:3
    %set(gca, 'YAxisLocation', 'origin') % plot results
    subplot(3,1,k)
    plot(1:(nIter+1), sun_error(:,k), '-');
    grid on;
    title('Error angle sun sensor vs sun direction with albedo');
    %legend(lab(k),strcat(lab(k), ' estimated'));
   
end
    end

