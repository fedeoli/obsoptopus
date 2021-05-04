%% Geometric method
% This method compute the position estimate for the entire  fleet. The
% computation is therefore centralized. 
% Input: 
%   1) Chi: (nagent x 3) matrix containing the a priori estimate for all
%   the agents.
%   2) GPS: (nagent x 3) matrix containing the GPS measures for all the
%   agents.
%   3) adjmat_UWB: (nagent x nagent) matrix ontaining the set of relative
%   distances measured by UWB. 
%   4) theta: first parameter
%   5) beta: second parameter
%   6) Check_dist: boolean flag used to trigger the security check on the
%   final estimate. 
% Output:
%   1) opt: data structure containing the set of nagent position estimates.
function opt = Position_opt_cloud_num_v5(Chi, GPS, adjmat_UWB, theta, beta, check_dist)
%     fprintf("Optimizing - geometric method\n");
    
    % nagent
    nagent = size(Chi,1);
    
    % safe check
    dist_thresh = 5e-2;
    
    % Versor Vij computation - GPS values: using the GPS values the direction
    % between agent i and agent j is computed. Data are stored in a 3D matrix.
    % Index (i,j,k) defines the kth component of the vector connecting agent i
    % to agent j.
    direction_matrix_Chi = zeros(nagent,nagent,3);
    for i = 1:nagent
       for j = 1:nagent
          if i ~= j
            direction_matrix_Chi(i,j,:) = (Chi(j,:) - Chi(i,:))/norm(Chi(j,:) - Chi(i,:));
          end
       end
    end 

    % Position according to UWB measures. Agent i position can be derived
    % from agent j assuming to know Dji and Vji. Thus, each agent's position
    % can be estimated by all the others, resulting in nagent-1 different
    % estimates. 
    Chi_UWB = zeros(nagent,nagent,3);
    for i = 1:nagent
        for j = 1:nagent
            if i ~= j
                for k = 1:3
                    Chi_UWB(i,j,k) = Chi(j,k) + direction_matrix_Chi(j,i,k)*adjmat_UWB(j,i);
                end
            end
        end
    end

    % New position estimate. All the estimated positions of agent i from agents
    % j are averaged, obtaining the centroid of the UWB position points cloud.
    % This value is averaged once more with agent i GPS position measured by
    % agent i itself. 
    Chi_estimate = zeros(nagent,3);
    centroid = zeros(nagent,3);
    for i = 1:nagent
        for k = 1:3
            centroid(i,k) = sum(Chi_UWB(i,:,k))/(nagent -1);
            Chi_estimate(i,k) = beta*Chi(i,k) + (1-beta)*(((1-theta)*centroid(i,k) + theta*GPS(i,k)));
        end  
    end
    
    % check if the new position is too far from the starting point
    if check_dist == 1
        for i = 1:nagent
           if  (norm(Chi(i,:)-Chi_estimate(i,:)) >= dist_thresh)
               Chi_estimate(i,:) = Chi(i,:);
           end
        end
    end
   
    opt.Chi_est = Chi_estimate;
    % optional
%     opt.centroid = centroid;
%     opt.Chi_UWB = Chi_UWB;
%     opt.direction_matrix_Chi = direction_matrix_Chi;
end


