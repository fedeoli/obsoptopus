%% Attitude_Obsever_V2_1: estimates the i-th satellite attitude using MAHONY
function  [xnew, Pnew, map_out] = Observer_EKF_v1(map, params,x,z)

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

%%% Handling a priori attitude coordinates
%%% Apriori Estimated Quaternion Dynamics as in Mahony 10.1109/TAC.2008.923738
[xhat, params] = map.model_propagate(map.ActualTimeIndex-1,map.Ts,x,params);

%%% estimated measure %%%
measure_forward = 1;
[map.buf_dyhat, Yhat] = map.get_measure(xhat,0,measure_forward,map.buf_dyhat,map.intYhat_full_story,params,map.ActualTimeIndex);
map.Yhat_full_story(:,end+1) = Yhat(:,1);
map.dYhat_full_story(:,end+1) = Yhat(:,2);
map.intYhat_full_story(:,end+1) = Yhat(:,3);
z_hat = Yhat(:,1);

% position offset
offset = map.integration_pos*6;

%%%%% EKF %%%%%
[x_hat_next, Pnext] = ekf_V6(map,params,xhat,z,z_hat);

% update cvariance matrix
Pnew = Pnext;

% update estimation state + measures
xhat(1:offset) = x(1:offset);
xhat(offset+1:offset+4) = x_hat_next(1:4);
xhat(offset+5:offset+7) = x_hat_next(5:7);

% estimation assignment
xnew = xhat;

map_out = map;



end

