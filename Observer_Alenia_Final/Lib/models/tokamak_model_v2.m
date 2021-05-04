function x_dot = tokamak_model_v2(t,x,params)

    % Dynamic matrix
    params.A = [0, 1,0,0; 0,0,1,0; ...
                0,-params.c0+params.gamma0,-params.c1,params.gamma1; ...
                -params.Ki*params.gainPHD*params.beta,-params.Kp*params.gainPHD*params.beta,-params.Kd*params.gainPHD*params.beta,-params.beta];
    x_dot = 5.*params.A*x;
end