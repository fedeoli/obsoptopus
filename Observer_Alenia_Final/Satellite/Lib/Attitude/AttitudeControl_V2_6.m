function tau = AttitudeControl_V2_6(satellites_attitude, satellites_iner_ECI, pos, params)

%   AttitudeControl_V2_4.m
%   Made by Sapienza Gn Lab
%
%   Computes the attitude control torque to be applied to each deputy to make it track (reach) the desired attitude.
%
%   INPUT
%   satellites_attitude: array (7*(N+1) x 1) containing the current attitude and angular velocities of each satellite (N is the number of deputy satellites)
%   satellites_iner_ECI: Array (6*(N+1) x 1) containing the satellites' (chief + deputies) inertial coordinates (N is the number of deputy
%                        satellites).
%   params: structure containing the following fields
%       - params.Ndeputy = number of deputy satellites
%       - params.control_dir_Hill = matrix (3 x N) containing the orbital control of each deputy
%       - params.ControlStrategy = flag for the selection of the Control Strategy to be adopted (either 'Impulsive' or 'Continuous')
%       - params.ParkingTime = time instant to which the control activation can be delayed
%       - params.DesiredAttitude = matrix (3 x N) containing the euler angles representing the desired attitude for each deputy based on the firing direction (N is the
%                    number of deputy satellites)
%       - params.mi = Earth's gravitational constant [km^3*s^-2]
%       - params.kp = Proportional attitude control gain
%       - params.kd = Derivative attitude control gain
%
%   OUTPUT
%   tau: matrix (3 x (N+1)) containing the attitude torque to be applied to each satellite (N is the number of deputy satellites). The torque of the
%        chief is always set to zero.
%
%   VERSION
%   20191125 V2_1:
%   - First Release
%
%   20200113 V2_2:
%   - With Impulsive control strategy: the attitude control is computed between two subsequent impulses to follow the required mission profile.
%
%   20200519 V2_3:
%   - Conversions between different attitude sets replaced by 'RotationConversion_V2_1', a purposely developed, non built-in Matlab function
%   - The attitude control is computed also if t < params.ParkingTime & params.CollisionFlag = 1, i.e. if a collision occurs before
%     the end of the parking time.
%
%   %   20200908 V2_4:
%   - Fixed bug for which the 'FiringAttitude' field of 'params.sat(i)' structure was used before it was defined when using the continuous control strategy. 
%     In particular, if the control did not fire because under the firing threshold, the firing attitude was not defined yet, but called in this function.

%%%%%%%%%%%%%
global DynOpt
%%%%%%%%%%%%%

% time
t = DynOpt.time(pos);

% Extract the number of deputy satellites
Ndeputy = params.Ndeputy;

% Initialize the output matrix
tau = zeros(3, Ndeputy + 1);

% Compute ECI2Hill rotation matrix and angular velocity
vect_r = satellites_iner_ECI(1:3);
vect_v = satellites_iner_ECI(4:6);
R_ECI2Hill = RECI2Hill(vect_r, vect_v);
omega_Hill2ECI_ECI = OmegaHill2ECI(vect_r, vect_v, params.mi);

% Extract the desired attitude for each deputy
DesiredAttitude = params.DesiredAttitude;

% Compute the attitude torque for each satellite
for i = 1:Ndeputy + 1
    
    if (t > params.ParkingTime && i > 1 && ( ( strcmpi(params.ControlStrategy, 'Continuous') && ( params.u_module(i) > 0 || ( isfield(params.sat(i), 'FiringAttitude') && isempty(params.sat(i).FiringAttitude) == 0 ) ) ) || ( strcmpi(params.ControlStrategy, 'Impulsive') && params.sat(i).CanManeuverAttitude == 1 ) ) ) || (i == 1 && params.sat(1).CanManeuverAttitude == 1) || ( t <= params.ParkingTime && params.CollisionFlag == 1 && params.sat(i).CA_Maneuvering == 1 )
         
        % Convert the desired attitude to dcm representation
        R_Hill2Body = RotationConversion_V2_1('EA321toDCM', DesiredAttitude(:,i)'*180/pi);
        
        % Refer the desired attitude wrt ECI reference frame and convert to quaternions
        R_ECI2Body = R_Hill2Body*R_ECI2Hill;
        
        % convert reference trajectory
        q_ref = RotationConversion_V2_1('EA321toQ', DesiredAttitude(:,i)'*180/pi);
        q_ref = ConvertQuat_V2_1(q_ref, 'ScalarTo1');
        
        % store quat ref
        DynOpt.quat_ref(:,t) = q_ref;
        
        % Extract the deputy current attitude and angular velocity
        q = satellites_attitude(1 + 7*(i-1): 4 + 7*(i-1))';
        
        % get omega 
        omega_Body2ECI_Body = satellites_attitude(5 + 7*(i-1): 7 + 7*(i-1));
        
        % testing
        omega_ref = R_ECI2Body*omega_Hill2ECI_ECI;
%         omega_ref = zeros(3,1);
        
        
        % Compute the deputy current angular velocity wrt to Hill reference frame
        omega_Body2Hill_Body = omega_Body2ECI_Body - omega_ref;
        
        % Compute the quaternion error between current and desired attitudes
        q_err = QuaternionError_V2_1(q, q_ref);
        
        % Compute PD attitude control torque
        tau(:,i) = -sign(q_err(1))*params.sat(i).kp*q_err(2:4)' - params.sat(i).kd*omega_Body2Hill_Body;
        
    end
    
end

end