%% EKF A matrix numeric %%
function A = Amatrix_EKF_v2(map,params,x,u)

    % quaternion
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    
    % omega
    wx = x(5);
    wy = x(6);
    wz = x(7);
    
    % magnetometers
    if map.nMagneto >= 1
        Bx1 = map.mag_field_vector(1);
        By1 = map.mag_field_vector(2);
        Bz1 = map.mag_field_vector(3);
    end
    if map.nMagneto >= 2
        Bx2 = map.mag_field_vector(1);
        By2 = map.mag_field_vector(2);
        Bz2 = map.mag_field_vector(3);
    end
    
    %%% params %%%
    Ixx = params.sat(1).I(1,1);
    Iyy = params.sat(1).I(2,2);
    Izz = params.sat(1).I(3,3);
    mie = params.mi;
    r = params.SatellitesCoordinates(1:3);
    rx = r(1);
    ry = r(2);
    rz = r(3);
    r_norm = norm(r);
    v = params.SatellitesCoordinates(1:3);
    vx = v(1);
    vy = v(2);
    vz = v(3);
    v_norm = norm(v);
    
    % input
    taux = u(1);
    tauy = u(2);
    tauz = u(3);
    
    % magnetometers
    dt1 = map.RPYbetweenMagSensors(1);
    dt2 = map.RPYbetweenMagSensors(2);
    dt3 = map.RPYbetweenMagSensors(3);
    
    
    %%% matrix computation %%%
    A = map.Asym(Ixx, Iyy, Izz, ...
                    wx, wy, wz, ...
                    q1, q2, q3, q0, ...
                    rx, ry, rz, r_norm, ...
                    vx, vy, vz, v_norm, ...
                    mie, ...
                    Bx1, By1, Bz1, ...
                    Bx2, By2, Bz2, ...
                    dt1, dt2, dt3, ...
                    taux, tauy, tauz);
                
    A = double(A);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
end