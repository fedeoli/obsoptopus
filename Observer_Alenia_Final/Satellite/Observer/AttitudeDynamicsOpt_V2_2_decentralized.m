 function [dw]  = AttitudeDynamicsOpt_V2_2_decentralized(w, params,ith_satellite,k_microstep)
%A global observer for attitude and gyro biases from vector measurements
% Philippe Martin  Ioannis Sarras
%10.1016/j.ifacol.2017.08.1868

%alpha hat = w(1:3)
%beta hat = w(4:6)
%b omega hat = w(7:9)

global ObserverTest Agent

% Initialize integration vectors and parameters
dw = zeros(length(w), 1);
alpha_hat = w(1:3);
beta_hat = w(4:6);
bias_omega_hat = w(7:9);

alpha_m = Agent(ith_satellite).magneto(:,ObserverTest.actual_time_index); %sampled
beta_m = Agent(ith_satellite).magneto2(:,ObserverTest.actual_time_index); %sampled
gyros_m = Agent(ith_satellite).gyros(:,ObserverTest.actual_time_index);

dw(1:3) = cross(alpha_hat,gyros_m-bias_omega_hat) - ObserverTest.attitude_Km*(alpha_hat-alpha_m);
dw(4:6) = cross(beta_hat,gyros_m-bias_omega_hat) - ObserverTest.attitude_Kmbis*(beta_hat-beta_m);
dw(7:9) = 0*(0.1*cross(alpha_hat,alpha_m) + 0.1*cross(beta_hat,beta_m));
