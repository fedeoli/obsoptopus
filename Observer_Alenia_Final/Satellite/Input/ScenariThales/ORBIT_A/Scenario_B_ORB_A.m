%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     ANALYSIS SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the final animation
animation_on = 0;                           % animation_on = 1 :  A final animation of the whole formation is displayed

% Activation of .gif file creation
createGIF = 0;                              % createGIF = 1 : A .gif file of the final animation is created. Only works if animation_on = 1

% Activation of the GPS noise
GPSnoise = 0;                               % GPSnoise = 1 : the deputy relative states used for control calculation is obtained by GPS measures

% Activation of the Control
params.Control = 0;                         % params.Control = 1 : The Control is applied to the formation
                                            % params.Control = 0 : No Control is applied to the formation, only free dynamics is integrated
                                            % IMPORTANT NOTE: The value of the Control flag does not influence the presence of the collision avoidance
                                            % strategy (see also 'params.CollisionAvoidance' flag comment)

% Selection of the Impulsive Control Strategy
params.ControlStrategy = 'Impulsive';       % params.ControlStrategy = 'Impulsive' : Selects the Impulsive Control Strategy
                                            % params.ControlStrategy = 'Continuous' : Selects the Continuous Control Strategy

% Activation of the Collision Probability calculation
params.CollisionProbabilityCalculation = 0; % params.CollisionProbabilityCalculation = 1 : The Collision Probability is computed at each iteration
                                            % params.CollisionProbabilityCalculation = 0 : The Collision Probability is not computed at any iteration
                                          
% Activation of the Collision Avoidance strategy
params.CollisionAvoidance = 0;              % params.CollisionAvoidance = 1 : The Collision Avoidance control calculation is activated, regardless of the value 
                                            %                                 set for the params.Control flag
                                            % params.CollisionAvoidance = 0 : The Collision Avoidance control calculation is deactivated

% Disable the Collison Avoidance strategy if the Collision Probability is not being calculated                                            
if ~params.CollisionProbabilityCalculation
    
    params.CollisionAvoidance = 0;
    
end

% Enable Collision Probability calculation if Collison Avoidance is selected
if params.CollisionAvoidance
    
    params.CollisionProbabilityCalculation = 1;
    
end

% Selection of the Collision Avoidance Strategy (works only if 'params.CollisionAvoidance = 1')
params.CollisionAvoidanceStrategy = 'SafetyParking';         % params.CollisionAvoidanceStrategy = 'Potential' : An exponential control is applied based on the repulsive potential between two 
                                                             %                                                   deputies. This strategy can only be adopted if 'Continuous' control strategy is selected. 
                                                             % params.CollisionAvoidanceStrategy = 'SafetyImpulse' : The colliding deputy exerts an optimal Delta V which guarantees a minimum 
                                                             %                                                       separation distance at the expected collision location. As soon
                                                             %                                                       as the Delta V is applied, the deputy can resume its nominal mode.
                                                             % params.CollisionAvoidanceStrategy = 'SafetyDrift' : The colliding deputy exerts an optimal Delta V which guarantees a minimum 
                                                             %                                                     separation distance at the expected collision location. After having applied such impulse,
                                                             %                                                     the colliding deputy keeps drifting on the new established orbit for a fixed amount of time
                                                             %                                                     (which can be set by means of the 'params.CA_Time' parameter) after which it
                                                             %                                                     can resume its nominal mode.
                                                             % params.CollisionAvoidanceStrategy = 'SafetyParking' : The colliding deputy exerts an optimal Delta V which sends him on a
                                                             %                                                       safe parking orbit for a fixed amount of time (which can be set by means of the 
                                                             %                                                       'params.CA_Time' parameter) after which it can resume its 
                                                             %                                                       nominal mode.
  
% Activation of the global CA safe mode (works only if 'params.CollisionAvoidance = 1')
params.CA_GlobalSafeMode = 1;               % params.CA_GlobalSafeMode = 1 : If a Collision is spotted ahead in time, each deputy stops maneuvering and enters CA_mode
                                            % params.CA_GlobalSafeMode = 0 : If a Collision is spotted ahead in time, the deputy which is asked to
                                            %                                perform the CA maneuver does it, while the others do not enter in CA_mode 
                                            %                                (e.g., they do not stop maneuvering if they were) 
                                            % IMPORTANT NOTE: The 'GlobalSafeMode' flag does not work with the 'Potential' CA strategy. 
                                            
% Activation of the two-axes-control (only for Continuous Control)
params.twoAxesControl = 0;                  % params.twoAxesControl = 1 : the deputy satellites are only able to apply the control force on y and z axis

% Activation of Augmented LQR
params.ALQR = 0;                            % params.ALQR = 1 : the control is evaluated by means of the Augmented LQR strategy (Riccati differential equation is solved)
                                            % params.ALQR = 0 : the control is evaluated by means of the LQR strategy (Riccati algebraic equation is solved)

% Activation of J2 model in ALQR control evaluation (only for Continuous Control - only if "params.ALQR = 1")
params.ALQR_J2 = 0;                         % params.ALQR_J2 = 1 : the Riccati differential equation is solved with the model taking into account J2 linearized effects
                                            % params.ALQR_J2 = 0 : the Riccati differential equation is solved with the model taking into account the
                                            %                      linearized dynamics without perturbations

% Activation of the FeedForward Control (only for Continuous Control)
params.FF = 0;                              % params.FF = 1 : the control is evaluated by means of both feedback and feedforward terms
                                            % params.FF = 0 : the control is evaluated by means of the sole feedback term
                                            % IMPORTANT NOTE: if Augmented LQR is selected, the Feedforward term is included in control evaluation,
                                            %                 regardless of the value assigned to "params.FF".
 
% Activation of the TV Observer
params.Observer = 0;                        % params.Observer = 1 : The Observer is activated
                                            % params.Observer = 0 : The Observer is not activated
                                           
% Activation of the attitude
params.Attitude = 0;                        % params.Attitude = 1 : The attitude dynamics, control and its influence on orbital maneuvering is included
                                            % params.Attitude = 0 : The attitude dynamics, control and its influence on orbital maneuvering is not included

% Activation of the real thrusters
params.RealThruster = 0;                    % params.RealThruster = 1 : A real thruster is included in the actuation of the control
                                            % params.RealThruster = 0 : The thrust is applicated "as is". No real thruster is being considered
                                            
                                            
% Format of the plot window
params.PlotWindowFormat = 'mixed';       % params.PlotWindowFormat = 'single' : plots in a single big window (two separate windows for orbital and atttiude plots)
                                            % params.PlotWindowFormat = 'multiple' : each quantity is plotted in a separate window

% Format of the time vector in the output
params.PlotTimeFormat = 'seconds';          % params.PlotTimeFormat = 'seconds' : time is plotted in seconds
                                            % params.PlotTimeFormat = 'hours' : time is plotted in hours
                                                                                              

% Selection of the active perturbation
J2_on = 1;                                  % J2_on = 1 :  J2 perturbation acting on the formation
Drag_on = 1;                                % Drag_on = 1 :  Drag perturbation acting on the formation
params.Drag_on = Drag_on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    IMPORTING CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Constants_V1_1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                CHIEF SATELLITE CHARACTERISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ChiefGeometry_V2_3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   CHIEF ORBIT DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ChiefOrbit_ORB_A;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 CHIEF ATTITUDE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ChiefAttitude_V2_2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SIMULATION TIMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_step = 1;                                                         % integration step [s]. IMPORTANT NOTE: If a high thrust engine is simulated, the time 
                                                                       %                       step should be decreased accordingly. With Busek BGT-X5, the
                                                                       %                       recommended value is time_step = 0.1 sec. With GomSpace CGP3 and 
                                                                       %                       Busek BIT-1, the recommended value is time_step = 1 sec.              

if params.Attitude == 0 && params.RealThruster == 0 
    
    time_step = 20;
    
end

tfinal = 3*86400;                                                         % simulation's final time [s]
tfinal = 500;
t = 0;                                                                  % time initialization
time = 0:time_step:tfinal;                                              % initialization of time array
tlength = length(time);
params.tfinal = tfinal;
params.time_step = time_step;

% Time to wait until control is applied to the chaser(s)
params.ParkingTime = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   DEPUTIES ORBITS DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the invariant condition
J2invariant = 0;                            % J2invariant = 1 :  if CI_MODE = 3 has been chosen, then the deputy's orbit is set to be J2 invariant


% Selection of the initial condition mode
CI_MODE = 1;                                % CI_MODE = 1 :  Initial conditions given in terms of relative coordinates: (x0, y0, z0, x0_dot , y0_dot , z0_dot)
                                            % CI_MODE = 2 :  Initial conditions given in terms of shape of the closed trajectory: (rho1, rho2, rho3, alpha0, beta0)
                                            % CI_MODE = 3 :  Initial conditions given in terms of relative orbital elements: (Da, De, Di, Dom, DOM, Dth)

% Deputy's inertial position and velocity calculation based on intial condition mdoe selction
switch CI_MODE
    
    case 1
        
        % Scenario D1 Specifications (Baseline cross-track 400m)
        x01 = 0.5;
        y01 = 1;
        z01 = 0.5;
        alpha01 = 0;
        x02 = 1;
        y02 = 2;
        z02 = 1;
        alpha02 = 120*pi/180;
        x03 = 1.5;
        y03 = 3;
        z03 = 1.5;
        alpha03 = 240*pi/180;
        
        % Initial conditions given in terms of relative coordinates: (x0, y0, z0, x0_dot , y0_dot , z0_dot)
        deputy_rel0_LVLH(1, 1:6) = [x01*cos(alpha01); -y01*sin(alpha01); z01*sin(alpha01); -x01*n*sin(alpha01); -y01*n*cos(alpha01); z01*n*cos(alpha01)];    % second deputy                                                                                                        % first deputy
        deputy_rel0_LVLH(2, 1:6) = [x02*cos(alpha02); -y02*sin(alpha02); z02*sin(alpha02); -x02*n*sin(alpha02); -y02*n*cos(alpha02); z02*n*cos(alpha02)];    % second deputy
        deputy_rel0_LVLH(3, 1:6) = [x03*cos(alpha03); -y03*sin(alpha03); z03*sin(alpha03); -x03*n*sin(alpha03); -y03*n*cos(alpha03); z03*n*cos(alpha03)];    % third deputy
        
        
        
        N_deputy = size(deputy_rel0_LVLH, 1);                           % number of deputy satellites
        chief_iner_ECI = satellites_iner_ECI;
        
        for i = 1:N_deputy
            
            % Transformation from relative to inertial coordinates
            satellites_iner_ECI(1 + 6*i : 6*(i + 1), 1)  = rel2iner_V2_2(deputy_rel0_LVLH(i, 1:6), chief_iner_ECI, chief_OOE, params);
            
        end
     
        
    case 2
        
        % Initial conditions given in terms of shape of the closed trajectory: (rho1, rho2, rho3, alpha0, beta0)
        deputy_shape(1, 1:5) = [1, 0, 0.5, pi/2, pi/2];                         % first deputy
        % deputy_shape(2, 1:5) = [1, 0, 2, 2*pi/3, 2*pi/3];                         % second deputy
        % deputy_shape(3, 1:5) = [1, 0, 2, 4*pi/3, 4*pi/3];                         % third deputy
        N_deputy = size(deputy_shape, 1);                                         % number of deputy satellites
        chief_iner_ECI = satellites_iner_ECI;
        
        for i = 1: N_deputy
            
            % Transformation from (rho1, rho2, rho3, alpha0, beta0) to inertial coordinates
            satellites_iner_ECI(1 + 6*i : 6*(i + 1), 1)  = shape2iner_V1_1(deputy_shape(i, 1:5), chief_iner_ECI, chief_OOE, params);
            
        end
        
        deputy_rel0_LVLH = AbsECI2RelHill_V1_2(satellites_iner_ECI, mi)';
        
        
    case 3
        
        % Initial conditions given in terms of relative_orbital elements: (Da, De, Di, Dom, DOM, Dth)
        deputy_dcoe(1,:) = [a_tandem - a_terrasar, e_tandem - e_terrasar, incl_tandem - incl_terrasar, argp_tandem - argp_terrasar, raan_tandem - raan_terrasar, 0];
        % deputy_dcoe(2,:) = [0, 0, 0.001, 0.001, 0, 0];                              % second deputy
        N_deputy = size(deputy_dcoe,1);                                             % number of deputy satellites
%         deputy_dcoe(1,4) = 0;
        deputy_dcoe(1,1) = 0;
%         deputy_dcoe(1,5) = 0;
        
        for i = 1:N_deputy
            
            if J2invariant
                
                deputy_MOE = J2invariantCOE_V1_1(chief_MOE, deputy_dcoe(i,:), params);
                deputy_OOE  = moe2ooe_V1_1(deputy_MOE, params);
                deputy_dcoe(i,:) = deputy_OOE - chief_OOE;
                
            end
            
            % Transformation from differential orbital elements to inertial coordinates
            satellites_iner_ECI(1 + 6*i : 6*(i + 1), 1) = dcoe2iner_V1_1(deputy_dcoe(i,:), chief_OOE, params);
            
        end
        
        deputy_rel0_LVLH = AbsECI2RelHill_V1_2(satellites_iner_ECI, mi)';
        
end

params.Ndeputy = N_deputy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   DEPUTIES GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DeputiesGeometry_V2_3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   DEPUTIES ATTITUDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DeputiesAttitude_V2_2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CONTROL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Engine %%%%%%%
GomSpaceCGP3;

%%%%%%% Attitude control parameters %%%%%%%
AttitudeControlParams_V2_1;

%%%%%%% Orbital control parameters %%%%%%%
OrbitalControlParams_V2_2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                COLLISION AVOIDANCE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CollisionAvoidanceParams_V2_2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                COMPLETE SATELLITE STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Build the 'sat' structure %%%%%%

% First Satellite (chief)
params.sat(1).CanManeuverAttitude = 1;
params.sat(1).MissionAttitude = 'Nadir pointing x';                 % Each satellite can perform attitude maneuvers which are aimed to accomplish a certain task.
                                                                    % The chief can always be set to maneuver its attitude, regardless of the scenario involved
                                                                    % and the flag set for the simulation (the only flag needed to activate the attitude maneuvering is 
                                                                    % params.Attitude = 1).
                                                                    % The deputy satellites, instead, can perform maneuver attitudes only between impulses, so only if
                                                                    % params.ControlStrategy = 'Impulsive'.
                                                                    % Possible mission attitude inputs:
                                                                    % - Nadir pointing x: The x Body axis points towards Earth
                                                                    % - Nadir pointing y: The y Body axis points towards Earth
                                                                    % - Nadir pointing z: The z Body axis points towards Earth 
                                                                    % - No maneuver: The attitude control is switched off (i.e., between impulses for deputies)
                                                                    % - Hold Attitude: The attitude is controlled to track the last firing attitude (this mode only works
                                                                    %   with deputies)

% % Second satellite (first deputy)
% params.sat(2).TrajectoryProfile = @HelixWithPhase_V1_1;
% params.sat(2).xRef = x0;
% params.sat(2).yRef = y0;
% params.sat(2).zRef = z0;
% params.sat(2).Phase = 0;
% params.sat(2).MissionAttitude = 'Nadir Pointing x';                 
% params.sat(2).IsWorking = 1;
% 
% % Third satellite (second deputy)
% params.sat(3).TrajectoryProfile = @HelixWithPhase_V1_1;
% params.sat(3).xRef = x0;
% params.sat(3).yRef = y0;
% params.sat(3).zRef = z0;
% params.sat(3).Phase = 120;
% params.sat(3).MissionAttitude = 'Nadir Pointing x';
% params.sat(3).IsWorking = 1;
% 
% % Fourth satellite (third deputy)
% params.sat(4).TrajectoryProfile = @HelixWithPhase_V1_1;
% params.sat(4).xRef = x0;
% params.sat(4).yRef = y0;
% params.sat(4).zRef = z0;
% params.sat(4).Phase = 240;
% params.sat(4).MissionAttitude = 'Nadir Pointing x';
% params.sat(4).IsWorking = 1;

% Fill the 'sat' structure
FillSatelliteStructure_V2_2;