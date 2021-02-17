%% Init section

% Symbolic definitions
syms p_x p_y p_z;        % Agent position
syms v_x v_y v_z;        % Agent velocity

% global structure containing EKF pos symbolic data
global EKF_pos

%Set params ObserverTest

% chief_OOE = [[6988.44162521818,0.00174269961850053,1.70815092407859,
% 1.62815615642515,4.92249794527405e-08,-0.0577608933031420]
chief_OOE = params.chief_OOE;

% symbolic state variables
myGPS_pos = [p_x p_y p_z];
myGPS_vel = [v_x v_y v_z];

% Vettore di STATO
X = [myGPS_pos, myGPS_vel];

% Vett. MISURA
h = [myGPS_pos];

% data structure assignment
EKF_pos.myGPS_pos = myGPS_pos;
EKF_pos.myGPS_vel = myGPS_vel;
EKF_pos.X = X;
EKF_pos.h = h;
EKF_pos.chief_OOE = chief_OOE;

