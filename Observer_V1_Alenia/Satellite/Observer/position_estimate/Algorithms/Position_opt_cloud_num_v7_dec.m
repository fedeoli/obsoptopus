%% Geometric method
% This method compute the position estimate for a single agent.
% Input: 
%   1) Chi: (nagent x 3) matrix containing the a priori estimate for all
%   the agents.
%   2) GPS: (1 x 3) matrix containing the GPS measures for the agent
%   considered.
%   3) adjmat_UWB: (nagent x nagent) matrix ontaining the set of relative
%   distances measured by UWB. 
%   4) j_select: number of the agent considered.
%   5) theta: first parameter
%   6) beta: second parameter
%   7) Check_dist: boolean flag used to trigger the security check on the
%   final estimate. 
% Output:
%   1) opt: data structure containing the set of nagent position estimates.
function opt = Position_opt_cloud_num_v7_dec(Chi, GPS, adjmat_UWB, j_select, theta, beta, check_dist, projection)
%     fprintf("Optimizing - geometric method\n");

    global Agent
    
    % nagent
    nagent = size(Chi,1);
    
    % safe check
    dist_thresh = 5e-2;
    
    % packet loss
    trust_agent_UWB = ones(1,nagent);
    
    % Versor Vij computation - Chi values: using the Chi values the direction
    % between agent i and agent j is computed. Data are stored in a matrix.
    % Index (i,k) defines the kth component of the vector connecting agent i
    % to agent j_select.
    direction_matrix_Chi = zeros(nagent,3);
    for i = 1:nagent
        if Agent(j_select).SuccessfullyReceivedData(i)
            if i ~= j_select 
                if strcmp(projection,'Chi')
                    direction_matrix_Chi(i,:) = (Chi(j_select,:) - Chi(i,:))/norm(Chi(j_select,:) - Chi(i,:));
                elseif strcmp(projection,'GPS')
                    direction_matrix_Chi(i,:) = (GPS(j_select,:) - GPS(i,:))/norm(GPS(j_select,:) - GPS(i,:));
                end
            end
        else
            trust_agent_UWB(i) = 0;
        end
    end 

    % Position according to UWB measures. Agent i position can be derived
    % from agent j assuming to know Dji and Vji. Thus, each agent's position
    % can be estimated by all the others, resulting in nagent-1 different
    % estimates.
    Chi_UWB = zeros(nagent,3);
    for i = 1:nagent
        if Agent(j_select).SuccessfullyReceivedData(i) 
            if i ~= j_select
                for k = 1:3
                    if strcmp(projection,'Chi')
                        Chi_UWB(i,k) = Chi(i,k) + direction_matrix_Chi(i,k)*adjmat_UWB(j_select,i);
                    elseif strcmp(projection,'GPS')
                        Chi_UWB(i,k) = GPS(i,k) + direction_matrix_Chi(i,k)*adjmat_UWB(j_select,i);
                    end
                end
            end
        end
    end

    % nagent trusted
    nagent_trusted = length(find(trust_agent_UWB == 1));
    if trust_agent_UWB(j_select) == 1
       nagent_trusted = nagent_trusted - 1; 
    end
    
    % Create trustful Chi matrix
    consider_centroid_flag = 0;
    Chi_UWB_select = zeros(nagent_trusted,3);
    for i = 1:nagent
        if Agent(j_select).SuccessfullyReceivedData(i)
            if i ~= j_select
                if i > j_select
                    Chi_UWB_select(i-1,:) = Chi_UWB(i,:);
                else
                    Chi_UWB_select(i,:) = Chi_UWB(i,:);
                end
                consider_centroid_flag = 1;
            end
        end
    end
    % remove zero columns
    Chi_UWB_select( all(~Chi_UWB_select,2), : ) = [];
    
    % New position estimate. All the estimated positions of agent j_select from agents
    % j are averaged, obtaining the centroid of the UWB position points cloud.
    % This value is averaged once more with agent i GPS position measured by
    % agent itself. 
    Chi_estimate = zeros(1,3);    
    centroid = zeros(1,3); 
    if consider_centroid_flag == 1
        for k = 1:3
            centroid(k) = sum(Chi_UWB_select(:,k))/nagent_trusted;
        end
    end
    
    % decide whether UWB shall be used or not
    if consider_centroid_flag == 1  %  && nagent_trusted > 1 (CHOOSE)
        UWB_flag = 1;
    else
        UWB_flag = 0;
    end
    
    % single agent position estimate
    for k = 1:3
        if UWB_flag
            if strcmp(projection,'Chi')
                Chi_estimate(k) = beta*Chi(j_select,k) + (1-beta)*(((1-theta)*centroid(k) + theta*GPS(k)));
            elseif strcmp(projection,'GPS')
                Chi_estimate(k) = beta*Chi(j_select,k) + (1-beta)*(((1-theta)*centroid(k) + theta*GPS(j_select,k)));
            end
        else
            if strcmp(projection,'Chi')
                Chi_estimate(k) = beta*Chi(j_select,k) + (1-beta)*GPS(k);
            elseif strcmp(projection,'GPS')
                Chi_estimate(k) = beta*Chi(j_select,k) + (1-beta)*GPS(j_select,k);
            end
        end
    end
    
    % check if the new position is too far from the starting point
    if check_dist == 1
       if  (norm(Chi(j_select,:) - Chi_estimate) >= dist_thresh)
               Chi_estimate = Chi(j_select,:);
       end
    end
   
    opt.Chi_est = Chi_estimate;
    % optional 
%     opt.centroid = centroid;
%     opt.Chi_UWB = Chi_UWB;
%     opt.direction_matrix_Chi = direction_matrix_Chi;
end


