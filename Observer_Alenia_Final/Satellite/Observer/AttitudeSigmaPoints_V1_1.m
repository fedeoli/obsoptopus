function [W0m,W0c,Wim,CHIi] = AttitudeSigmaPoints_V1_1(x_m,P_x)

%  Daniele Carnevale (31.03.2020), implementation of 
%  Vinther, K., Fuglsang Jensen, K., Larsen, J. A., & Wisniewski, R. (2011). Inexpensive CubeSat Attitude
%  Estimation Using Quaternions and Unscented Kalman Filtering. Automatic Control in Aerospace, 4(1).
% Differences: here the first component of the quaternion is the fourth in
% the paper.
%  x_m: previous estimated state [quaternion, omegas, bias_mag, bias_gyro],  the first component of the quaternion is the scalar part
%  P_x: previous covariance matrices
%--------------------------------
%testing parameters
% b_mag = 0.1;
% b_gyro = 0.1;
% P_x = blkdiag(eye(3),eye(3),b_mag*eye(3),b_gyro*eye(3))*1E-3;
% x_m = [[1 1 0 0] [2 1 -1 ] [2 1 1] [1 2 -3] ]; 
%--------------------------------


L = numel(x_m)-1; %the quaternion has only three degree of freedom that's way -1, size(P) = (L-1)x(L-1)
alpha = sqrt(3); %as in the paper
beta = 2; %as in the paper
cappa = 1; %as in the paper
lambda = 1;%alpha^2*(L+cappa)-L; %scaling parameter, as in the paper

CHIi = zeros(L+1,2*L+1);
deltaCHIi = zeros(L, 2*L);
q_prev = x_m(1:4)'; %previous quaternion
omegas_prev = x_m(5:7)';
mag = x_m(8:10)';
gyro =  x_m(11:13)';

%Weights
W0m = lambda/(L+lambda);
W0c = lambda/(L+lambda) + (1-alpha^2+beta);
Wim = 1/(2*(L+lambda));

try
    A = chol(P_x);
catch
    disp('Attitude Sigma Point function: issue with Choleskii factorization of P')
    P_x = (P_x + P_x')/2.0; %1E-3*eye(size(P_x))
    A = chol(P_x);
    %eig(A)
end

temp = (sqrt(L+lambda)*A')';
deltaCHIi = [ zeros(L,1), temp, -temp ]; %1.1: error signa points, the
%uncorrect selection that might yield imaginary quaternions

%TO SOLVE THE ISSUE OF NON UNIT QUATERNION
%paper A Simple Attitude Unscented Kalman Filter: Theory and Evaluation in a Magnetometer-Only Spacecraft Scenario
%MURTY S. CHALLA et al., Digital Object Identifier 10.1109/ACCESS.2016.2559445
%DeltaTheta = temp(:,1:3);
%DeltaOmega = temp(:,4:6);
%DeltaMag = temp(:,7:9);
%DeltaGyro = temp(:,10:12);

for ind=1:2*L+1,
    %if((deltaCHIi(1:3,ind)'*deltaCHIi(1:3,ind))>1)
    %    disp('ERROR: AttitudeSigmaPoints_V1_1, ||deltaCHIi(1:3,ind)||>1');
    %end
    %CHIi(:,ind) = [ (quatmultiply([sqrt(1-(deltaCHIi(1:3,ind)'*deltaCHIi(1:3,ind))) deltaCHIi(1:3,ind)' ], q_prev))...
     %                    (omegas_prev + deltaCHIi(4:6,ind)'  )  (mag + deltaCHIi(7:9,ind)')  (gyro + deltaCHIi(10:12,ind)') ]';
    
     CHIi(1:4,ind) = quatnormalize(quatmultiply( q_prev, [ 1, deltaCHIi(1:3,ind)'/2.0 ] ));
     CHIi(5:7,ind) =   (omegas_prev' + deltaCHIi(4:6,ind) );
     CHIi(8:10,ind) =   (mag' + deltaCHIi(7:9,ind) );
     CHIi(11:13,ind) =   (gyro' + deltaCHIi(10:12,ind) );
end


