%% Attitude_Obsever_V2_1: estimates the i-th satellite attitude using MAHONY
function  [xnew, Pnew, map_out] = Observer_EKF_v3(map, params, x_past, x_now, z_now, z_hat, u_now)

%%%% INIT SECTION %%%%
%%%% extract if position is integrated %%%%
% position offset
offset = map.integration_pos*6;

% only attitude
xhat_now = x_now(offset+1:offset+7,end);

%%%% Linearisation %%%%
% Linearized State equation in xk-1
G = Gmatrix_EKF_v2(map,params,x_past,u_now); 
% Linearized State equation in xk 
H = Hmatrix_EKF_v2(map,params,x_now,u_now,z_now);

%%%% reset covariance %%%%
if mod(map.ActualTimeIndex,map.Preset) == 0
    map.P = map.P0;
end

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
xnew = x_now;
xnew(offset+1:offset+4) = quatnormalize(transpose(x_hat_new(1:4)));
xnew(offset+5:offset+7) = x_hat_new(5:7);

% output map 
map_out = map;



end

