function EKF_position_V1(params,i,tspan_att,k)
%% EKF to estimate position and speeds of the agents

% Author: Daniele Giorgi, Federico Oliva
% Date: 15/09/2020

% global vars
global Agent ObserverTest EKF_pos

% symbolic vars
myGPS_pos = EKF_pos.myGPS_pos;
myGPS_vel = EKF_pos.myGPS_vel;
X = EKF_pos.X;
h = EKF_pos.h;
chief_OOE = EKF_pos.chief_OOE;


%RESET OF P after a given amount of samples
if(mod(i+1,ObserverTest.position_P_reset_aftersamples)==0)
    %disp(['P Reset - det(P)=  '   num2str(det(Agent(k).P)) ' substituted with det(Preset) = ' num2str(det(reshape(ObserverTest.PiAll(Ntest,:),size(ObserverTest.Pi))))])
    Agent(k).P = reshape(ObserverTest.PiAll(ObserverTest.Ntest,:),size(ObserverTest.Pi));
end

% symbolic covariance matrix
Agent(k).P_pos = Agent(k).P;

%% INITIAL CONDITIONS

% GET MEASURES from GPS 
if(Agent(k).UWBcorrectionActive(i+1)==1)
    myGPS = Agent(k).GPSopt(:,i+1);
else
    myGPS = Agent(k).GPS(:,i+1);
end

% MEASUREMENT assigment
%ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS,j) = prime 3 colonne
%del j-esimo agent, ossia le posizioni del GPS
% for j=1:2*ObserverTest.Na+1
%     if(ObserverTest.AddDerivativeEstimation ~= 0)% && i >ObserverTest.PseudoDerivative_d)
%         ZETA(:,j) = [ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS,j);ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS+3,j)];
%     else
%         ZETA(:,j) = [ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS,j)];
%     end
% end

%% Dynamic section

%Mat linearizzata STATO (con J2 effect)
A = AJ2_V1_1(i, chief_OOE, params); 

% compute output mapping jacobian
H = jacobian(h,X);

%observed state
% xhat_past = Agent(k).iner_ECI(1:3,i)';
xhat_past = Agent(k).xHatUKF(:,i)';
xhat =  Agent(k).xHatUKF(:,i+1)';

% numeric evaluation of the output mapping
h_num = subs(h, myGPS_pos, xhat(1:3));
%This time the H linearization matrix is always constant, as a blkdiag of
%ones and zeros!

% measurement assignment - eq. 5.12 Report_EKF_global
z = myGPS';
% estimation - eq. 5.13 Report_EKF_global
z_hat = double(h_num);

% Params storage procedure
params.P = Agent(k).P_pos;
params.Q = Agent(k).Q_pos;
params.R = Agent(k).R_pos;
params.R_ECI2Body = Agent(k).R_ECI2Body;
params.X = X;
params.A = A;
params.H = H;
params.ObserverTest = ObserverTest;
params.tspan_att = tspan_att;

%%%%% EKF %%%%%
[xHat, Pnext] = ekf_V5(xhat_past,xhat,z,z_hat,params,0,0);

% post filter assignments
if(i < ObserverTest.OptimizationsEachSampling ) %it is not the last that need to update everything
    Agent(k).xHatUKFmultiple = xHat;
else
    % update cvariance matrix
    Agent(k).P = Pnext;
    Agent(k).xHat = xHat; %the first component is the variations with respect to r
    Agent(k).xHatUKF(:,i+1) = Agent(k).xHat;
    Agent(k).xaUKF = [Agent(k).xHat; zeros(ObserverTest.Ndisturbance,1)]; %[xHat; zeros(Ndisturbance,1)]; %note that xaUKF contains the estimated state at time i+1-Sensor_delay
    
    Agent(k).iner_ECI_EstimationError(:,i+1) = Agent(k).iner_ECI(:,i+1) - Agent(k).xHatUKF(:,i+1);
    Agent(k).NormEstimationErrorXYZ(i+1) = norm(Agent(k).iner_ECI_EstimationError(1:3,i+1));
    ObserverTest.estimated_satellites_iner_ECI(1+(k-1)*6:6+(k-1)*6,i+1) = Agent(k).xHatUKF(:,i+1);
    
end
end
