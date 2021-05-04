%% EKF filter method
% EKF   Extended Kalman Filter for nonlinear dynamic systems
% [x, P] = ekf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P
% for nonlinear dynamic system:
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R

% INPUT: 
% 1) state vector - previous time instant
% 2) state vector - current time instant
% 3) measurements vector - current time instant
% 4) estimated measurement vector - current time instant
% 5) parameters
% 6) number of magnetometers used in the system (max 2)
% 7) number of Sun sensors used in the system (max 1)
% OUTPUT
% 1) a posteriori estimated state vector
% 2) a posteriori estimated covariance matrix
function [x_hat_next, P]= ekf_V5(xk0,xk1,yk,z_hat,params,nMagneto,Sun)

% global vars
global ObserverTest
%omega_Body2ECI_Body = params.omega_Body2ECI_Body;
%q_ECI2Body = params.q_ECI2Body;
%attitude_ECI = params.attitude_ECI;
P = params.P;
Q = params.Q;
R = params.R;
X = params.X;
A = params.A;
H = params.H;
BI_sym1 = params.BI1_sym;
if nMagneto == 2
   BI_sym2 = params.BI2_sym;
end
if Sun == 1
   PSun = params.Psun; 
end
R_ECI2Body = params.R_ECI2Body;
tspan = params.tspan_att;

% %Set params ObserverTest
% ObserverTest = params.ObserverTest;

% Linearized State equation in xk-1
A0 = double(subs(A, X, xk0)); 

% Convert magnetometer measures to inertial frame (quaternion dependence)
if nMagneto == 1
    MagTemp_I1 = R_ECI2Body'*yk(1:3)';
elseif nMagneto == 2
    MagTemp_I1 = R_ECI2Body'*yk(1:3)';
    MagTemp_I2 = R_ECI2Body'*yk(4:6)';
end
if Sun == 1
    h_meas = R_ECI2Body'*yk(end-5:end-3)';
end

% Linearized State equation in xk
H1 = subs(H, X, xk1);
if nMagneto == 1
    H1 = subs(H1, BI_sym1, MagTemp_I1');
elseif nMagneto == 2
    H1 = subs(H1, BI_sym1, MagTemp_I1');
    H1 = subs(H1, BI_sym2, MagTemp_I2');
end
if Sun == 1
    H1 = subs(H1, PSun, h_meas');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a priori covariance - solve riccati equation
% P0 = P;
% [~,temp] = ode45(@(t,P)mRiccati(t,P,A0,Q), tspan, P0);
% 
% % reshape output
% temp = reshape(temp(end,:),size(A));
% P0 = temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = size(A0,1);
%%%%%%%%%%%%% A priori covariance - simple method %%%%%%%%%%%%%%%%%%%%%%
phi = eye(dim)+A0*(tspan(2)-tspan(1));
P0 = phi*P*phi'+ Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kalman gain
K = P0*H1'*(pinv(double(H1*P0*H1' + R)));

amp = 1e0;
K = double(amp*K);
delta_attitude_xHatMeno = K*(yk - z_hat)';

% CHECK PRINTS
% size(K)
% double(P)   
% double(A0)
% double(P0)
% double(H1')
% pinv(double(H1*P0*H1' + R))
% double(P0*H1')
% norm(double(K)) 
% err = (yk'- z_hat')
% delta_attitude_xHatMeno

% CONTROLLA CHE I QUATERNIONI SIANO SEMPRE NORMALIZZATI AD OGNI STEP

%state estimate CORRECTION
x_hat_next= xk1' + delta_attitude_xHatMeno; 

% AVOIDING OF NON-NORMALIZED QUATERNIONS (Challa)
if (ObserverTest.EKF_attitude)
    q_aposteriori = quatnormalize(x_hat_next(1:4)');
    x_hat_next(1:4) = q_aposteriori;
end

% covariance estimation
temp = (eye(dim) - K*H1)*P0;

%NEW state covariance matrix
P = double(temp);      
end