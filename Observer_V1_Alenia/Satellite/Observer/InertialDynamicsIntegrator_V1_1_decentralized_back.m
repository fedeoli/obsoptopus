function dx = InertialDynamicsIntegrator_V1_1_decentralized_back(stato, params,ith_satellite,k_microstep)

% back in time for the observer 
%everything for the satellite i (i-th)


global  fh_c r_c h_c fd_c ObserverTest

% Extraction of constants from "params" structure
mi = params.mi;             % Earth's planetary constant
Re = params.Re;             % Earth's equatorial mean radius
J2 = params.J2;             % Earth's J2 coefficient
omega_e = params.omega_e;   % Earth's rotational velocity
CDAF = params.CDAF;         % Satellite's balistic coefficients
u_LVLH = reshape(ObserverTest.params_u_alltimes(:,ObserverTest.u_time_index),3,(ObserverTest.Nagents-1));          % Control matrix

% Number of satellites
n_sat = 1;
% Redefinition of control matrix, including a zero control for the chief satellite
u_LVLH = [zeros(3,1), u_LVLH];

% Arrays initialization
dx = zeros(6*n_sat, 1);
fn_ECI = zeros(3, n_sat);
p = zeros(n_sat, 1);
r = zeros(n_sat, 1);
th = zeros(n_sat, 1);
f2r = zeros(n_sat, 1);
f2th = zeros(n_sat, 1);
f2h = zeros(n_sat, 1);
fj2_ECI = zeros(3, n_sat);
pos = zeros(3, n_sat);
vel = zeros(3, n_sat);
v_rel = zeros(3, n_sat);
v_rel_mod = zeros(n_sat, 1);
rho = zeros(n_sat, 1);
D = zeros(n_sat, 1);
fd_ECI = zeros(3, n_sat);
u_ECI = zeros(3, 1);

% Transformation form inertial position and velocity to classical orbital elements
coe(1:6) = rv2coe_V1_1(stato(1:3), stato(4:6), mi);

% Newton's gravitational force computation
fn_ECI = LVLH2ECI_V1_1([-mi/norm(stato(1:3))^2,0,0], coe(3), coe(5), coe(4) + coe(6));

% J2 perturbation computation
p = coe(1)*(1-coe(2)^2);                                                                                     % semilatum rectum of the i-th satellite
r = p/(1 + coe(2)*cos(coe(6)));                                                                           % radius of the i-th satellite
th = coe(6) + coe(4);                                                                                        % argument of the latitude of the i-th satellite
f2r = -3/2*J2*mi*Re^2/(r^4)*(1-(3*(sin(coe(3)))^2*(sin(th))^2));                     % radial component of J2 perturbation acting on i-th satellite
f2th = -3/2*J2*mi*Re^2/(r^4)*(sin(coe(3)))^2*sin(2*th);          % tangential component of J2 perturbation acting on i-th satellite
f2h = -3*J2*mi*Re^2/(r^4)*sin(coe(3))*cos(coe(3))*sin(th);   % out-of-plane component of J2 perturbation acting on i-th satellite
fj2_ECI= LVLH2ECI_V1_1([f2r, f2th, f2h], coe(3), coe(5), coe(4) + coe(6))';  % J2 perturbation acting on i-th satellite (expressed in ECI)

% Drag perturbation computation
pos = stato(1 : 3);
vel = stato( 4 : 6);
v_rel = [vel(1) + pos(2)*omega_e; vel(2) - pos(1)*omega_e; vel(3)];                                              % satellite's relative velocity wrt Earth's rotation
v_rel_mod = norm(v_rel);                                                                                                % module of the relative velocity
%if(r - Re < 0)
%    disp('InertialDynamicsIntegrator_V1_1_decentralized: r - Re < 0');
%    Re
%    r
%end
%r-Re
rho = ExponentialAtmDensity_V1_1(r - Re);                                                                               % atmospheric density at the altitude of i-th satellite
D = -(1/2)*CDAF(ith_satellite)*rho*v_rel_mod^2;                                                                                    % drag acceleration module
fd_ECI(1) = D*v_rel(1)/v_rel_mod;
fd_ECI(2) = D*v_rel(2)/v_rel_mod;
fd_ECI(3) = D*v_rel(3)/v_rel_mod;

% Control Reconstructed as the true one, i.e. the one given to the target plant 
chief_coe_control = rv2coe_V1_1(ObserverTest.StateChief4Control(1:3,k_microstep), ObserverTest.StateChief4Control(4:6,k_microstep), mi);
u_ECI = LVLH2ECI_V1_1(u_LVLH(:,ith_satellite), chief_coe_control(3), chief_coe_control(5), chief_coe_control(4) + chief_coe_control(6));   % Rotation of control vector from LVLH (chief) to ECI

%%DD
% State derivatives computation  - Absolute dynamics ECI
dx(1) = -stato(4);
dx(2) = -stato(5);
dx(3) = -stato(6);
dx(4) = -(fn_ECI(1) + fd_ECI(1) + fj2_ECI(1) + u_ECI(1)*ObserverTest.ControlOn);
dx(5) = -(fn_ECI(2) + fd_ECI(2) + fj2_ECI(2) + u_ECI(2)*ObserverTest.ControlOn);
dx(6) = -(fn_ECI(3) + fd_ECI(3) + fj2_ECI(3) + u_ECI(3)*ObserverTest.ControlOn);

%?
%dx = 0*dx;


% if( ith_satellite == 1) %only for the estimated chief satellite
%     chief_coe = coe;
%     r_c = p(1)/(1 + chief_coe(2)*cos(chief_coe(6)));
%     h_c = sqrt(mi*p(1));
%     fd_c = fd_ECI(:,1);
%     TotalPerturbation_ECI = [(fj2_ECI(1) + fd_c(1)), (fj2_ECI(2) + fd_c(2)), (fj2_ECI(3) + fd_c(3))];                             % Sum of perturbation expressed in ECI reference frame
%     TotalPerturbation_LVLH = ECI2LVLH_V1_1(TotalPerturbation_ECI, chief_coe(3), chief_coe(5), (chief_coe(4) + chief_coe(6)));                                   % Sum of perturbations expressed in LVLH reference frame
%     fh_c = TotalPerturbation_LVLH(3);
% end
% 
end