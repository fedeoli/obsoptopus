function dx = InertialDynamicsIntegrator_V2_1_decentralized(stato, params,ith_satellite,k_microstep)

%   InertialDynamicsIntegrator_V2_1.m

global  ObserverTest Agent

% Extraction of constants from "params" structure
mi = params.mi;             % Earth's planetary constant
Re = params.Re;             % Earth's equatorial mean radius
J2 = params.J2;             % Earth's J2 coefficient
omega_e = params.omega_e;   % Earth's rotational velocity
u_LVLH = params.u;          % Control matrix

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
fd_Body = zeros(3, n_sat);
fd_ECI = zeros(3, n_sat);
u_ECI = zeros(3, n_sat);

% Loop over the number of satellites
%for i = 1:n_sat
    
    % Transformation form inertial position and velocity to classical orbital elements
    coe(1:6) = rv2coe_V1_1(stato(1:3), stato(4:6), mi);

    % Newton's gravitational force computation
    fn_ECI(:) = LVLH2ECI_V1_1([-mi/norm(stato(1:3))^2,0,0], coe(3), coe(5), coe(4) + coe(6));

    % J2 perturbation computation
    p = coe(1)*(1-coe(2)^2);                                                                                     % semilatum rectum of the i-th satellite
    r = p/(1 + coe(2)*cos(coe(6)));                                                                           % radius of the i-th satellite
    th = coe(6) + coe(4);                                                                                        % argument of the latitude of the i-th satellite
    f2r = -3/2*J2*mi*Re^2/(r^4)*(1-(3*(sin(coe(3)))^2*(sin(th))^2));                                               % radial component of J2 perturbation acting on i-th satellite
    f2th = -3/2*J2*mi*Re^2/(r^4)*(sin(coe(3)))^2*sin(2*th);                                                        % tangential component of J2 perturbation acting on i-th satellite
    f2h = -3*J2*mi*Re^2/(r^4)*sin(coe(3))*cos(coe(3))*sin(th);                                             % out-of-plane component of J2 perturbation acting on i-th satellite
    fj2_ECI(:) = LVLH2ECI_V1_1([f2r, f2th, f2h], coe(3), coe(5), coe(4) + coe(6))';      % J2 perturbation acting on i-th satellite (expressed in ECI)

    % Drag perturbation computation
    q_ECI2Body = Agent(ith_satellite).attitude_xHatUKF(1:4,ObserverTest.u_time_index);                                                  % Extract i-th satellite's attitude quaternions
    R_ECI2Body = quat2dcm(q_ECI2Body');                                                                                              % Compute i-th satellite's attitude matrix
    pos(:) = stato(1:3);
    vel(:) = stato(4:6);
    v_rel(:) = [vel(1) + pos(2)*omega_e; vel(2) - pos(1)*omega_e; vel(3)];                                              % satellite's relative velocity wrt Earth's rotation
    v_rel_body = R_ECI2Body*v_rel(:);                                                                                                  % Compute velocity in BRF
    v_rel_mod = norm(v_rel(:));                                                                                                % module of the relative velocity
    rho = ExponentialAtmDensity_V1_1(r - Re);                                                                                 % atmospheric density at the altitude of i-th satellite
    
     for j = 1:size(params.sat(ith_satellite).aero_prop, 2)
        
        % Extract i-th satellite j-th surface properties
        Area = params.sat(ith_satellite).aero_prop(j).A;
        Surf_normal = params.sat(ith_satellite).aero_prop(j).n';
        cross_section = max(Area*dot(v_rel_body, Surf_normal), 0);                                                                  % j-th surface cross section
        fd_Body(:) = fd_Body(:) - (1/2)*(params.sat(ith_satellite).CD*cross_section/params.sat(ith_satellite).M)*rho*v_rel_mod^2*v_rel_body/v_rel_mod;  % Drag acceleration on j-th surface
        
     end
         
    fd_ECI(:) = R_ECI2Body'*fd_Body(:);                                         
    
    % Control Computation
    u_ECI(:) = LVLH2ECI_V1_1(u_LVLH(:,ith_satellite), ObserverTest.EstimatedChiefCoe(3), ObserverTest.EstimatedChiefCoe(5), ObserverTest.EstimatedChiefCoe(4) + ObserverTest.EstimatedChiefCoe(6));                                                       % Rotation of control vector from LVLH (chief) to ECI
    
    % State derivatives computation
    dx(1) = stato(4);
    dx(2) = stato(5);
    dx(3) = stato(6);
    dx(4) = fn_ECI(1) + fd_ECI(1) + fj2_ECI(1) + u_ECI(1);
    dx(5) = fn_ECI(2) + fd_ECI(2) + fj2_ECI(2) + u_ECI(2);
    dx(6) = fn_ECI(3) + fd_ECI(3) + fj2_ECI(3) + u_ECI(3);

%end

if(ith_satellite == 1)
   
    ObserverTest.r_c = p/(1 + coe(2)*cos(coe(6)));
    ObserverTest.h_c = sqrt(mi*p(1));
    ObserverTest.fd_c = fd_ECI(:);
    TotalPerturbation_ECI = [(fj2_ECI(1) + ObserverTest.fd_c(1)), (fj2_ECI(2) + ObserverTest.fd_c(2)), (fj2_ECI(3) + ObserverTest.fd_c(3))];                             % Sum of perturbation expressed in ECI reference frame
    TotalPerturbation_LVLH = ECI2LVLH_V1_1(TotalPerturbation_ECI, coe(3), coe(5), (coe(4) + coe(6)));                                   % Sum of perturbations expressed in LVLH reference frame
    ObserverTest.fh_c = TotalPerturbation_LVLH(3);
    ObserverTest.EstimatedChiefCoe = coe;

end