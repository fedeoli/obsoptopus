 function [dw]  = AttitudeDynamicsOpt_V1_1_decentralized(w, params,ith_satellite,k_microstep)
%Attitude integration function for attitude UKF, note that there are extra
%states as we try to estimate the bias of magnetometer and gyros as well.

global ObserverTest Agent

% Initialize integration vectors and parameters
dw = zeros(length(w), 1);
mie = params.mi;

% Compute Hill's attitude and angular velocity wrt ECI
%vect_r_chief = params.X_SAT_ECI(1:3);%Agent(1).xHatUKF(1:3,ObserverTest.actual_time_index);????
%vect_v_chief = params.X_SAT_ECI(4:6);%Agent(1).xHatUKF(4:6,ObserverTest.actual_time_index);????
vect_r_chief = Agent(1).xHatUKF(1:3,ObserverTest.actual_time_index);
vect_v_chief =Agent(1).xHatUKF(4:6,ObserverTest.actual_time_index);
R_ECI2Hill = RECI2Hill(vect_r_chief, vect_v_chief);
omega_Hill2ECI_ECI = OmegaHill2ECI(vect_r_chief, vect_v_chief, mie);

% Extract attitude, angular velocity and inertia of satellite "i"
Att_sat_i = w;
q_ECI2Body =  Att_sat_i(1:4);
omega_Body2ECI_Body = Att_sat_i(5:7);
I = params.sat(ith_satellite).I;


% Extract "i"-th satellite's inertial position and velocity
%vect_r = params.X_SAT_ECI(6*(ith_satellite-1) + 1 : 6*(ith_satellite-1) + 3);%Agent(ith_satellite).xHatUKF(1:3,ObserverTest.actual_time_index); %USE THE ESTIMATES ????
%vect_v = params.X_SAT_ECI(6*(ith_satellite-1) + 4 : 6*(ith_satellite-1) + 6);%Agent(ith_satellite).xHatUKF(4:6,ObserverTest.actual_time_index);
vect_r =  Agent(ith_satellite).xHatUKF(1:3,ObserverTest.actual_time_index); 
vect_v = Agent(ith_satellite).xHatUKF(4:6,ObserverTest.actual_time_index);

r = norm(vect_r);

%     omega_Body2Hill_Body = -(R_ECI2Body*omega_Hill2ECI_ECI - omega_Body2ECI_Body);

%%%%%%%%%%%%%%%%%% Attitude Kinematics Equations (Quaternions) %%%%%%%%%%%%%%%%%%
Om = [0, -omega_Body2ECI_Body(1), -omega_Body2ECI_Body(2), -omega_Body2ECI_Body(3);
    omega_Body2ECI_Body(1), 0, omega_Body2ECI_Body(3), -omega_Body2ECI_Body(2);
    omega_Body2ECI_Body(2), -omega_Body2ECI_Body(3), 0, omega_Body2ECI_Body(1);
    omega_Body2ECI_Body(3), omega_Body2ECI_Body(2), -omega_Body2ECI_Body(1), 0];
dw(1:4,1) = 0.5*Om*q_ECI2Body;


%%%%%%%%%%%%%%%%%%%%%%%% Attitude Dynamics Equations %%%%%%%%%%%%%%%%%%%%%%%%

R_ECI2Body = quat2dcm(q_ECI2Body');
Iom = I*omega_Body2ECI_Body;
cross_omIom = [omega_Body2ECI_Body(2)*Iom(3) - omega_Body2ECI_Body(3)*Iom(2);...
    omega_Body2ECI_Body(3)*Iom(1) - omega_Body2ECI_Body(1)*Iom(3);...
    omega_Body2ECI_Body(1)*Iom(2) - omega_Body2ECI_Body(2)*Iom(1)];

% Gravity Gradient Torque Computation
vers_o_Body = -R_ECI2Body*(vect_r/r);
Io = I*vers_o_Body;
cross_oIo = [vers_o_Body(2)*Io(3) - vers_o_Body(3)*Io(2);...
    vers_o_Body(3)*Io(1) - vers_o_Body(1)*Io(3);...
    vers_o_Body(1)*Io(2) - vers_o_Body(2)*Io(1)];
GG_torque = (3*mie/(r^3))*cross_oIo;

% Aerodynamic Drag Torque Computation
v_rel = [vect_v(1) + vect_r(2)*params.omega_e; vect_v(2) - vect_r(1)*params.omega_e; vect_v(3)];                        % satellite's relative velocity wrt Earth's rotation
v_rel_mod = norm(v_rel);                                                                                                % module of the relative velocity
rho = ExponentialAtmDensity_V1_1(norm(vect_r) - params.Re);                                                             % atmospheric density at the altitude of i-th satellite
CD = params.sat(ith_satellite).CD;                                                                                                  % extract i-th satellite CD coefficient
air_drag_torque = [0; 0; 0];                                                                                            % Initialize drag torque for i-th satellite

for j = 1:size(params.sat(ith_satellite).aero_prop, 2)
    
    % Extract i-th satellite j-th surface properties
    Area = params.sat(ith_satellite).aero_prop(j).A;
    Surf_normal = params.sat(ith_satellite).aero_prop(j).n;
    arm = params.sat(ith_satellite).aero_prop(j).arm;
    
    v_rel_body = R_ECI2Body*v_rel;                                                                                      % Compute velocity in BRF
    cross_section = max(Area*dot(v_rel_body, Surf_normal), 0);                                                          % j-th surface cross section
    F_drag = -(1/2)*rho*CD*v_rel_mod^2*cross_section*v_rel_body/v_rel_mod;                                              % Drag force on j-th surface
    T_drag = cross(arm, F_drag);                                                                                        % Drag torque caused by j-th surface
    air_drag_torque = air_drag_torque + T_drag;                                                                         % Sum j-th torque to already calculated ones
    
end

%tau = params.tau(:,ith_satellite);
tau = reshape(ObserverTest.tau(:,ObserverTest.u_time_index),3,ObserverTest.Nagents);
dw(5:7, 1) = I\( - cross_omIom + GG_torque + air_drag_torque + tau(:,ith_satellite));

%dynamics of the extra states

dw(8:10, 1) = zeros(3,1); %constant magnetometer bias
dw(11:13, 1) = zeros(3,1); %constant gyro bias (temperature compensated sensor)

