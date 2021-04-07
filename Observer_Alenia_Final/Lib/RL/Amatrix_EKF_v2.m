%% EKF A matrix numeric %%
function A = Amatrix_EKF_v2(map,params,x,u)

    % position offset
    offset = map.integration_pos*6;

    % quaternion
    q0 = x(offset+1);
    q1 = x(offset+2);
    q2 = x(offset+3);
    q3 = x(offset+4);
    
    % omega
    wx = x(offset+5);
    wy = x(offset+6);
    wz = x(offset+7);
    
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
    v = params.SatellitesCoordinates(4:6);
    vx = v(1);
    vy = v(2);
    vz = v(3);
    v_norm = norm(v);
    
    % input
    taux = u(1);
    tauy = u(2);
    tauz = u(3);
    
    
    %%% matrix computation %%%
    A = map.Asym(Ixx, Iyy, Izz, ...
                    wx, wy, wz, ...
                    q1, q2, q3, q0, ...
                    rx, ry, rz, r_norm, ...
                    vx, vy, vz, v_norm, ...
                    mie, ...
                    taux, tauy, tauz);
                
    A = double(A);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
end