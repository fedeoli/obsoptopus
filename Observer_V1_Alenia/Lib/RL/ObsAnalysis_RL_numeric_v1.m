%%%% NUMERIC ORC COMPUTATION %%%%
function dtheta = ObsAnalysis_RL_numeric_2M_v1(map,x,y)

    % quaternion
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    
    % magnetometers
    if map.nMagneto >= 1
        Bx1 = y(4);
        By1 = y(5);
        Bz1 = y(6);
    end
    if map.nMagneto >= 2
        Bx2 = y(7);
        By2 = y(8);
        Bz2 = y(9);
    end
    

    dtheta = [  0,   0,   0, 2.0*Bx1*q0 + 2.0*By1*q3 - 2.0*Bz1*q2, 2.0*By1*q0 - 2.0*Bx1*q3 + 2.0*Bz1*q1, 2.0*Bx1*q2 - 2.0*By1*q1 + 2.0*Bz1*q0, 2.0*Bx2*q0 + 2.0*By2*q3 - 2.0*Bz2*q2, 2.0*By2*q0 - 2.0*Bx2*q3 + 2.0*Bz2*q1, 2.0*Bx2*q2 - 2.0*By2*q1 + 2.0*Bz2*q0;
                0,   0,   0, 2.0*Bx1*q1 + 2.0*By1*q2 + 2.0*Bz1*q3, 2.0*Bx1*q2 - 2.0*By1*q1 + 2.0*Bz1*q0, 2.0*Bx1*q3 - 2.0*By1*q0 - 2.0*Bz1*q1, 2.0*Bx2*q1 + 2.0*By2*q2 + 2.0*Bz2*q3, 2.0*Bx2*q2 - 2.0*By2*q1 + 2.0*Bz2*q0, 2.0*Bx2*q3 - 2.0*By2*q0 - 2.0*Bz2*q1;
                0,   0,   0, 2.0*By1*q1 - 2.0*Bx1*q2 - 2.0*Bz1*q0, 2.0*Bx1*q1 + 2.0*By1*q2 + 2.0*Bz1*q3, 2.0*Bx1*q0 + 2.0*By1*q3 - 2.0*Bz1*q2, 2.0*By2*q1 - 2.0*Bx2*q2 - 2.0*Bz2*q0, 2.0*Bx2*q1 + 2.0*By2*q2 + 2.0*Bz2*q3, 2.0*Bx2*q0 + 2.0*By2*q3 - 2.0*Bz2*q2;
                0,   0,   0, 2.0*By1*q0 - 2.0*Bx1*q3 + 2.0*Bz1*q1, 2.0*Bz1*q2 - 2.0*By1*q3 - 2.0*Bx1*q0, 2.0*Bx1*q1 + 2.0*By1*q2 + 2.0*Bz1*q3, 2.0*By2*q0 - 2.0*Bx2*q3 + 2.0*Bz2*q1, 2.0*Bz2*q2 - 2.0*By2*q3 - 2.0*Bx2*q0, 2.0*Bx2*q1 + 2.0*By2*q2 + 2.0*Bz2*q3;
              1.0,   0,   0,                                    0,                                    0,                                    0,                                    0,                                    0,                                    0;
                0, 1.0,   0,                                    0,                                    0,                                    0,                                    0,                                    0,                                    0;
                0,   0, 1.0,                                    0,                                    0,                                    0,                                    0,                                    0,                                    0];
 
    
    
end