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
function [x_hat_next, P]= ekf_V6(map,params,xk,yk,z_hat)

% position offset
offset = map.integration_pos*6;

% get quaternion evolution
q_hat = xk(offset+1:offset+4,end)'; 
% get angular velocity evolution
omega_hat = xk(offset+5:offset+7,end)';
 
% state and oberved state
xhat = [q_hat, omega_hat]';

% Linearized State equation in xk-1
A = Amatrix_EKF_v2(map,params,xhat,yk); 

% Convert magnetometer measures to inertial frame (quaternion dependence)

% Linearized State equation in xk
H = Hmatrix_EKF_v2(map,params,xhat,yk);

%%%%%%%%%%%%% A priori covariance - simple method %%%%%%%%%%%%%%%%%%%%%%
dim = size(A,1);
P0 = A*map.P*A'+ map.Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kalman gain
K = P0*H'*(pinv(H*P0*H' + map.R));
K = double(K);

% innovation
delta_attitude_xHatMeno = K*(yk - z_hat);

% state estimate CORRECTION
x_hat_next= xhat + delta_attitude_xHatMeno; 

% covariance estimation
P = (eye(dim) - K*H)*P0;      
end