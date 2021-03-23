%% Attitude_Obsever_V2_1: estimates the i-th satellite attitude using MAHONY
function  [xnew, Pnew, map_out] = Observer_EKF_v2(map, params,x_past,z_now,u_now)

%%%% INIT SECTION %%%%
% position offset
offset = map.integration_pos*6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% propagate current state - S2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xhat_now_tot, params] = map.model_propagate(map.ActualTimeIndex-1,map.Ts,x_past,params);

%%% estimated measure %%%
measure_forward = 1;
[map.buf_dyhat, Yhat] = map.get_measure(xhat_now_tot,0,measure_forward,map.buf_dyhat,map.intYhat_full_story,params,map.ActualTimeIndex);
map.Yhat_full_story(:,end+1) = Yhat(:,1);
map.dYhat_full_story(:,end+1) = Yhat(:,2);
map.intYhat_full_story(:,end+1) = Yhat(:,3);
z_hat = Yhat(:,1);

%%%% extract if position is integrated %%%%
% get quaternion evolution
q_hat = xhat_now_tot(offset+1:offset+4,end)'; 
% get angular velocity evolution
omega_hat = xhat_now_tot(offset+5:offset+7,end)';
% state and oberved state
xhat_now = [q_hat, omega_hat]';

%%%% Linearisation %%%%
% Linearized State equation in xk-1
G = Amatrix_EKF_v2(map,params,x_past,u_now); 
% Linearized State equation in xk (u is useless here)
H = Hmatrix_EKF_v2(map,params,xhat_now,u_now);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% A priori covariance - S3 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pbar = G*map.P*G'+ map.Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Kalman gain - S4 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = Pbar*H'*(pinv(H*Pbar*H' + map.R));
K = double(K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% state estimate - S5 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_hat_new = xhat_now + K*(z_now - z_hat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% covariance estimation - S6 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pnew = (eye(size(G,1)) - K*H)*Pbar;    

% update estimation state + measures
xnew = xhat_now_tot;
xnew(offset+1:offset+4) = x_hat_new(1:4);
xnew(offset+5:offset+7) = x_hat_new(5:7);

% output map 
map_out = map;



end

