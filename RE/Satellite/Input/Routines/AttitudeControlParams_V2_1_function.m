function params = AttitudeControlParams_V2_1_function(params)
    % Attitude sensor
    params.AttitudeSensorNoise = 0.007*pi/180; % std deviation pointing accuracy
    params.MaxTorque = 0.004*1e-6; % max torque 0.004 Nm = 4e-9 kg*km/s2*km

    % Attitude actuation
    params.kp_deputy = 1e-9; % 6e-10; %3e-10;                               % Attitude control proportional gain (it is set equal for each deputy)
    params.kd_deputy = 1e-8; %2e-9;    % 4e-9;                             % Attitude control derivative gain (it is set equal for each deputy)
    params.kp_chief = 6e-5;
    params.kd_chief = 2e-4;
    params.att_threshold = 2;
end