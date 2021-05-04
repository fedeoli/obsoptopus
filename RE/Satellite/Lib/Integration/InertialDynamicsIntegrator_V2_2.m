function dx = InertialDynamicsIntegrator_V2_2(satellites_iner_ECI, params)

%   InertialDynamicsIntegrator_V2_2.m
%   Made by Sapienza Gn Lab
%
%   Integrates the inertial dynamics equation of N satellites. It includes Drag and J2 perturbations calculations. Control is applied to each deputy, i.e., not 
%   on the first satellite, which is the chief. This function has to be used inside an already built integration method (e.g. ode45 or a home-made rk4 
%   fixed-step method).
%
%   INPUT
%   satellites_iner_ECI: Array (6*N x 1), where N is the number of satellites. It contains satellites' inertial position and velocity.
%   params: Structure containing the following fields
%       - params.mi: Earth's planetary constant
%       - params.Re: Earth's equatorial mean radius
%       - params.J2: Earth's J2 coefficient
%       - params.omega_e: Earth's rotational velocity
%       - params.CDAF: satellites' balistic coefficients
%       - params.u: control matrix (3 x N). j-th column contains the control vector applied to j-th deputy
%       - params.SatellitesAttitude: matrix (7 x N) containing the attitude of each satellite at time "t".
%
%   OUTPUT
%   dx: Array (6*N x 1) containing the derivatives of the state. 
%
%   VERSION
%   20181219 V1_1:
%   - First Release
%
%   20191114 V1_2:
%   - Addition of the attitude to determine the firing direction
%
%   20200109 V2_1:
%   - Attitude of each satellite considered for the purpose of correct cross-section area calculation (drag perturbation)
%
%   20200521 V2_2:
%   - The drag acceleration computation is now independant from the satellites' attitude parameters if the params.Attitude flag is set to zero

% Extraction of constants from "params" structure
mi = params.mi;             % Earth's planetary constant
Re = params.Re;             % Earth's equatorial mean radius
J2 = params.J2;             % Earth's J2 coefficient
omega_e = params.omega_e;   % Earth's rotational velocity
u_LVLH = params.u;          % Control matrix

% Number of satellites
n_sat = length(satellites_iner_ECI)/6;

% Redefinition of control matrix, including a zero control for the chief satellite
% u_LVLH = [zeros(3,1), u_LVLH];
u_LVLH = [u_LVLH, u_LVLH];

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
for i = 1:n_sat
    
    % Transformation form inertial position and velocity to classical orbital elements
    coe(6*(i-1)+1:i*6) = rv2coe_V1_1(satellites_iner_ECI(6*(i-1)+1:6*(i-1)+3), satellites_iner_ECI(6*(i-1)+4:6*(i-1)+6), mi);

    % Newton's gravitational force computation
    fn_ECI(:,i) = LVLH2ECI_V1_1([-mi/norm(satellites_iner_ECI( 6*(i-1) + 1 : 6*(i-1) + 3 ))^2,0,0], coe(6*(i-1) + 3), coe(6*(i-1) + 5), coe(6*(i-1) + 4) + coe(6*(i-1) + 6));

    % J2 perturbation computation
    p(i) = coe(6*(i-1)+1)*(1-coe(6*(i-1)+2)^2);                                                                                     % semilatum rectum of the i-th satellite
    r(i) = p(i)/(1 + coe(6*(i-1)+2)*cos(coe(6*(i-1)+6)));                                                                           % radius of the i-th satellite
    th(i) = coe(6*(i-1)+6) + coe(6*(i-1)+4);                                                                                        % argument of the latitude of the i-th satellite
    f2r(i) = -3/2*J2*mi*Re^2/(r(i)^4)*(1-(3*(sin(coe(6*(i-1)+3)))^2*(sin(th(i)))^2));                                               % radial component of J2 perturbation acting on i-th satellite
    f2th(i) = -3/2*J2*mi*Re^2/(r(i)^4)*(sin(coe(6*(i-1)+3)))^2*sin(2*th(i));                                                        % tangential component of J2 perturbation acting on i-th satellite
    f2h(i) = -3*J2*mi*Re^2/(r(i)^4)*sin(coe(6*(i-1)+3))*cos(coe(6*(i-1)+3))*sin(th(i));                                             % out-of-plane component of J2 perturbation acting on i-th satellite
    fj2_ECI(:,i) = LVLH2ECI_V1_1([f2r(i), f2th(i), f2h(i)], coe(6*(i-1)+3), coe(6*(i-1)+5), coe(6*(i-1)+4) + coe(6*(i-1)+6))';      % J2 perturbation acting on i-th satellite (expressed in ECI)

    % Drag perturbation computation                                                                                       
    pos(:,i) = satellites_iner_ECI(6*(i-1) + 1 : 6*(i-1) + 3);
    vel(:,i) = satellites_iner_ECI(6*(i-1) + 4 : 6*(i-1) + 6);
    v_rel(:,i) = [vel(1,i) + pos(2,i)*omega_e; vel(2,i) - pos(1,i)*omega_e; vel(3,i)];                                              % satellite's relative velocity wrt Earth's rotation
    v_rel_mod(i) = norm(v_rel(:,i));                                                                                                % module of the relative velocity
    rho(i) = ExponentialAtmDensity_V1_1(r(i) - Re);                                                                                 % atmospheric density at the altitude of i-th satellite
    
    if params.Attitude
        
        q_ECI2Body = params.SatellitesAttitude(3*(i-1) + 1 : 3*(i-1) + 4);                                                              % Extract i-th satellite's attitude quaternions
        R_ECI2Body = quat2dcm(q_ECI2Body');                                                                                             % Compute i-th satellite's attitude matrix
        v_rel_body = R_ECI2Body*v_rel(:,i);                                                                                             % Compute velocity in BRF
        
        for j = 1:size(params.sat(i).aero_prop, 2)
            
            % Extract i-th satellite j-th surface properties
            Area = params.sat(i).aero_prop(j).A;
            Surf_normal = params.sat(i).aero_prop(j).n';
            cross_section = max(Area*dot(v_rel_body, Surf_normal), 0);                                                                           % j-th surface cross section
            fd_Body(:,i) = fd_Body(:,i) - (1/2)*(params.sat(i).CD*cross_section/params.sat(i).M)*rho(i)*v_rel_mod(i)^2*v_rel_body/v_rel_mod(i);  % Drag acceleration on j-th surface
            
        end
        
        fd_ECI(:,i) = R_ECI2Body'*fd_Body(:,i);
        
    else
                
        % Compute Drag acceleration in ECI reference frame
        fd_ECI(:,i) = -(1/2)*(params.sat(i).CD*params.sat(i).MeanCrossSection/params.sat(i).M)*rho(i)*v_rel_mod(i)*v_rel(:,i);
        
    end
    
    % Control Computation
    u_ECI(:,i) = LVLH2ECI_V1_1(u_LVLH(:,i), coe(3), coe(5), coe(4) + coe(6));                                                       % Rotation of control vector from LVLH (chief) to ECI
    
    % State derivatives computation
    dx(6*(i-1) + 1) = satellites_iner_ECI(6*(i-1) + 4);
    dx(6*(i-1) + 2) = satellites_iner_ECI(6*(i-1) + 5);
    dx(6*(i-1) + 3) = satellites_iner_ECI(6*(i-1) + 6);
    dx(6*(i-1) + 4) = fn_ECI(1,i) + fd_ECI(1,i) + fj2_ECI(1,i) + u_ECI(1,i);
    dx(6*(i-1) + 5) = fn_ECI(2,i) + fd_ECI(2,i) + fj2_ECI(2,i) + u_ECI(2,i);
    dx(6*(i-1) + 6) = fn_ECI(3,i) + fd_ECI(3,i) + fj2_ECI(3,i) + u_ECI(3,i);
    
    %%% SPEED UP MODIFICATION %%%
    dx = params.eps_coef*dx;

end


params.r_c = p(1)/(1 + coe(2)*cos(coe(6)));
params.h_c = sqrt(mi*p(1));
params.fd_c = fd_ECI(:,1);
TotalPerturbation_ECI = [(fj2_ECI(1,1) + params.fd_c(1)), (fj2_ECI(2,1) + params.fd_c(2)), (fj2_ECI(3,1) + params.fd_c(3))];                             % Sum of perturbation expressed in ECI reference frame
TotalPerturbation_LVLH = ECI2LVLH_V1_1(TotalPerturbation_ECI, coe(3), coe(5), (coe(4) + coe(6)));                                   % Sum of perturbations expressed in LVLH reference frame
params.fh_c = TotalPerturbation_LVLH(3);


end