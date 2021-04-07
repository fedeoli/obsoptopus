%% observability analysis on the nonlinear system (attitude)
function [theta,dtheta,dtheta_num] = ObsAnalysis_RL_v1(map,n_der,x,y,calc)

    if calc == 1
        % init
        X = map.X;
        f = map.f;
        h = map.h;
        dimState = length(X);
        dimOut = length(h);

        %%%%%%% COMPUTE THETA %%%%%%%
        n_iter = n_der;
        for k=1:n_iter
            if k==1
               % init theta
               theta = h; 
            else
               for z=1:dimOut
                   for i=1:dimState
                        theta_temp(z,i) = diff(theta(z+(k-2)*dimOut),X(i));
                   end
               end
               theta = vpa([theta; theta_temp*f],2); 
            end
        end

        %%%%%% COMPUTE DTHETA %%%%%%%%
        dimTheta = length(theta);
        for k=1:dimTheta
            for i=1:dimState
                dtheta(i,k) = diff(theta(k),X(i));
            end
        end
        dtheta = vpa(dtheta,2);
        dtheta_num = 0;
    else
        %%%%% COMPUTE RANK %%%%%%%%%
        dtheta_num = subs(map.dtheta,map.X,x);
        if map.nMagneto >= 1
            dtheta_num = subs(dtheta_num,map.magneto1,y(4:6));
        end
        if map.nMagneto >= 2
            dtheta_num = subs(dtheta_num,map.magneto2,y(7:9));
        end
        dtheta_num = double(dtheta_num);
        theta = 0;
        dtheta = 0;
    end
    
    
end

