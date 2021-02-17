%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   CHIEF ORBIT A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SENTINEL_1A
% Type: Near-Polar Sun-Synchronous
%
% Mean Local Solar Time at ascending node: 18:00 h (nominal)
% Repeat Cycle: 12 days
% Cycle Length: 175 orbits
% Inclination: 98.12 [deg]
% Altitude: 693 [km]
%
% TLE (celestrack)
% 1 39634U 14016A   20201.80753859 -.00000116  00000-0 -14878-4 0  9999
% 2 39634  98.1829 208.3314 0001272  85.7508 274.3829 14.59196703335281
%          incl    RAAN     ecc      argp    meanAnom revperday


% Chief's Mean Orbital Elements (MOE)
ecc = 0.0001272;                                                        % chief's orbit's mean eccentricity
i = 98.1829*pi/180;                                                     % chief's orbit's mean inclination [rad]
om = 85.7508*pi/180;                                                    % chief's orbit's mean argument of the perigee [rad]
RAAN = 208.3314*pi/180;                                                 % chief's orbit's mean right ascension of the ascending node [rad]
f_0 = 0*pi/180;                                                         % chief's mean true anomaly at time 0 [rad]
T = 86400/14.59196703335281;
n = 2*pi/T;                                                             % chief's orbit's mean motion [rad/s]
a = (mi/(n^2))^(1/3);                                                   % chief's orbit's mean semimajor axis [km]

chief_MOE = [a, ecc, i, om, RAAN, f_0];                                 % array with chief's mean orbital elements [A, E, I, ARG, RAAN, TH]
params.n = n;

% Transformation of chief's Mean Orbital Elements (MOE) to Osculating Orbital Elements (OOE)

chief_OOE = moe2ooe_V1_1(chief_MOE, params);
params.chief_OOE = chief_OOE;

% Conversion from chief's Osculating Orbital Elements to inertial coordinates expressed in ECI reference frame and storing into satellites_iner_ECI
satellites_iner_ECI = coe2rv_V1_1(chief_OOE(1:6), mi);