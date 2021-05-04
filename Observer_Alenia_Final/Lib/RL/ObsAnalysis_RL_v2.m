%% observability analysis on the nonlinear system (attitude)
function [theta,dtheta,dtheta_num] = ObsAnalysis_RL_v2(map,n_der,x,y,calc)

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
        if map.nMagneto >= 2
            dtheta_num = ObsAnalysis_RL_numeric_2M_v1(map,x,y);
        else
            dtheta_num = ObsAnalysis_RL_numeric_1M_v1(map,x,y);
        end
        theta = 0;
        dtheta = 0;
    end
    
    
end

