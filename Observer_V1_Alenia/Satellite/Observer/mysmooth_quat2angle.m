function [r1 r2 r3] = mysmooth_quat2angle( q,old_angles )
% Smoothing the  QUAT2ANGLE function
% old_angles are the past angle values
[r1,r2,r3] = quat2angle(q);

   %AVOIDING JUMPS
    if( abs(r1 - old_angles(1)) > 0.99*pi )
            r1 = r1 - sign(r1 - old_angles(1))*2*pi;
    end
    if( abs(r2 - old_angles(2)) > 0.99*pi )
            r2 = r2 - sign(r2 - old_angles(2))*2*pi;
    end
    if( abs(r3 - old_angles(3)) > 0.99*pi )
            r3 = r3 - sign(r3 - old_angles(3))*2*pi;
    end
    
end
