function [satellites_iner_ECI_estimated,...
               deputy_rel_LVLH_estimated,...
               satellites_attitude_estimated] = Observer_V2_5(satellites_iner_ECI, deputy_rel_LVLH, satellites_attitude, params,i,tspan)

%%  POSITION and ATTITUDE OBSERVER
%INPUT: 
% satellites_iner: Array (6*N x 1), with N equal to the number of satellites. It is built in the following way: [pos1; vel1; pos2; vel2; pos3; vel3;
%                    .... ], where the position and velocities are must be expressed in ECI reference frame. The first six components (pos1; vel1) are 
%                    the inertial coordinates of the chief satellite. 
% satellites_attitude: Array (7*N x 1), with N equal to the number of satellites. It is built in the following way: [quaternion1; omega1; quaternion2; omega2; etc...
%                    .... ], where the first component of the quaternion is    the scalar component
% deputy_rel_LVLH: deputy rLVLH coordinates with respect to chief
% params:  object defined by La Sapienza containing general costants,  simulation parameters, control torques, desired references...
% i: previous time step index (the estimates are computed at time t = (i+1)*time_step), IT IS ASSUMED THAT time_step_att = time_step  (position).
% tspan: Time interval used in the integration (prediction step) = [i*time_step, (i+1)*time_step]

%OTUPUT: 
%this function processes the global (objects) variables Agent and ObserverTest, updating their fields and returning the vectors of the estimated quantities for all agents, e.g. satellites_iner_ECI_estimated is a vector with 6*Nagents components (Nagents = 1 (chief) + Ndeputy) in the ECI reference stacked as: [[x y z dx dy dz]1-th agent,
%[x y z dx dy dz]2-th agent?. , [x y z dx dy dz] Nagents -th agent]
%The vector deputy_rel_LVLH_estimated contains LVLH coordinates of the deputies (Nagents-1), and the vector satellites_attitude_estimated contains [[quaternion]1-th agent, 
%[omega]1-th agent, [quaternion]2-th agent, [omega]2-th agent,  [quaternion] Nagents -th agent, [omega] Nagents -th agent].




global ObserverTest Agent 



%if the observer does not run then the true values are passed back
satellites_iner_ECI_estimated = satellites_iner_ECI;
deputy_rel_LVLH_estimated = deputy_rel_LVLH;
satellites_attitude_estimated = satellites_attitude;


ObserverTest.satellites_iner_ECI_alltimes(:,i+1) = satellites_iner_ECI;
ObserverTest.params_u_alltimes(:,i) = reshape(params.u,1,3*(ObserverTest.Nagents-1))';
ObserverTest.u_time_index = i;
ObserverTest.actual_time_index = i+1; 
%The control input are stored adding also the control 0 for the Chief
ObserverTest.u(:,i) = [zeros(3,1);reshape(params.u,3*(ObserverTest.Nagents-1),1)]; 
ObserverTest.tau(:,i) = reshape(params.tau,3*(ObserverTest.Nagents),1); %tau does contain already the attitude control for the Chief
    
    

%updating data not linked with estimation
for k=1:ObserverTest.Nagents,
    Agent(k).iner_ECI(:,i+1) = satellites_iner_ECI(1 + (k-1)*6:6 + (k-1)*6);
    Agent(k).attitude(:,i+1) = satellites_attitude(1 + (k-1)*7:7 + (k-1)*7);
    if(ObserverTest.ObserverON==0)
        Agent(k).xHatUKF(:,i+1) = Agent(k).iner_ECI(:,i+1);
    end
    if(ObserverTest.AttitudeObserverOn==0)
        Agent(k).attitude_xHatUKF(1:7,i+1) = Agent(k).attitude(:,i+1);
        Agent(k).attitude_xHatUKF(8:end,i+1) = [ObserverTest.MagnetoBias;ObserverTest.GyroBias];
    end
end



if(ObserverTest.ObserverON)
    
    ObserverTest.StateChief4Control = satellites_iner_ECI(1:6,:); %we store the state evolution within the integration step of the chief necessary to evaluate the control as the one really actuated
    
    if(ObserverTest.ThrusterEstimateUpdate == 1)
        %%%%%%%%%%%%%%%%%%%%%
        % Update the state in case thusters have been used in the previous time step
        
        % Compute ECI2Hill rotation matrix
        vect_r = Agent(1).xHatUKF(1:3,i);%satellites_iner_ECI(1:3);
        vect_v = Agent(1).xHatUKF(4:6,i);%satellites_iner_ECI(4:6);
        % Compute chief COE
        chief_coe = rv2coe_V1_1(satellites_iner_ECI(1:3), satellites_iner_ECI(4:6), params.mi);
        
        %for all deputies
        temp = size(params.DV(:, :, 1));
        if(ObserverTest.DVsizeMemo(2) < temp(2) )
            ObserverTest.DVsizeMemo =  temp;
            for k = 2:ObserverTest.Nagents,
                DV = params.DV(:, end, k-1);
                if(norm(DV)>0)
                    quat_ECI = Agent(k).attitude_xHatUKF(1:4,i)'; %satellites_attitude(1 + 7*j: 4 + 7*j)'; %the first one is the chief among my agents
                    R_ECI2Body = quat2dcm(quat_ECI);
                    R_ECI2Hill = RECI2Hill(vect_r, vect_v);
                    R_Hill2Body = R_ECI2Body*R_ECI2Hill';
                    
                    DV_applied = R_Hill2Body'*[norm(DV); 0; 0];
                    
                    deputy_rel_LVLH(:,k-1) = deputy_rel_LVLH(:,k-1) + [zeros(3,1); DV_applied];
                    
                    %Modifying the estimated j-th deputy inertial coordinates
                    Agent(k).xHatUKF(:,i) = rel2iner_V2_2(deputy_rel_LVLH(:,k-1)', satellites_iner_ECI(1:6)', chief_coe, params);
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%
    
    %Evaluation of the A PRIORI self estimate for each agent using the
    %dynamic model: note that we temporarly set the estimate at time
    %xHatUKF(:,i+1) TO DEFINE THE EXHANGED DATA and which is going to
    %be replaced by the proper UKF algorithm next
    APrioriEstimationXYZ = zeros(1,3*ObserverTest.Nagents);
    
%%%% COMMENTO FEDEOLI %%%%
UKF_flag = 1;
GPSopt_flag = 1;

if UKF_flag == 1
    for k = 1:ObserverTest.Nagents,
        Agent(k).iner_ECI(:,i+1) = satellites_iner_ECI((k-1)*6+1:(k-1)*6+6);
        X = rk4_V1_1_decentralized(@InertialDynamicsIntegrator_V2_1_decentralized, tspan, Agent(k).xHatUKF(:,i), params,k);
        Agent(k).xHatUKF(:,i+1) =  X(:,end);%satellites_iner_ECI([1:6]+6*(k-1)); %??????
        ObserverTest.APrioriEstimationXYZ(1+3*(k-1):3+3*(k-1)) = Agent(k).xHatUKF(1:3,i+1);
    end
else
    for k = 1:ObserverTest.Nagents,
        Agent(k).iner_ECI(:,i+1) = satellites_iner_ECI((k-1)*6+1:(k-1)*6+6);
        X = rk4_V1_1_decentralized(@InertialDynamicsIntegrator_V2_1_decentralized, tspan, Agent(k).iner_ECI(:,i), params,k);
        Agent(k).xHatUKF(:,i+1) =  X(:,end);%satellites_iner_ECI([1:6]+6*(k-1)); %??????
        ObserverTest.APrioriEstimationXYZ(1+3*(k-1):3+3*(k-1)) = Agent(k).xHatUKF(1:3,i+1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %MEASUREMENTS
    Observer_Measurements_V1_3(satellites_iner_ECI,deputy_rel_LVLH, params,i);
    
    
    %possible multiple observer iteration between samples
    iteration = 1;
    while (iteration <= ObserverTest.OptimizationsEachSampling)
        
        if GPSopt_flag == 1
            GPS_Optimization_V1_geometric; %GPS optimization
        end
        
        if UKF_flag == 1
            Position_UKF_V1_3; %UKF estimator for the satellites position
        end
        
        
        %distance error after the a posteriori estimation step
        ObserverTest.AposterioriAdjancyMatrix =  setAdjacencyMatrixNorm(ObserverTest.estimated_satellites_iner_ECI(ObserverTest.PositionArrayIndex,i+1),ObserverTest.Nagents);
        for k=1:ObserverTest.Nagents,
            Agent(k).MismatchAposterioriDistances(i+1) = sum(abs(ObserverTest.AposterioriAdjancyMatrix(k,:)-ObserverTest.MeasuredDistances(k,:) ));
            Agent(k).ErrorAposterioriDistances(i+1) = sum(abs(ObserverTest.AposterioriAdjancyMatrix(k,:)-ObserverTest.AdjacencyMatrix(k,:) ));
        end
        
        iteration = iteration + 1;
        if(ObserverTest.UWBoptimizationOn==0) %stop if no UWB optimization is selected, there is no meaning in performing multiple GPS optimization!
            iteration = ObserverTest.OptimizationsEachSampling + 1;
        end
    end
    
    % CONTROL EVALUATION FOR THE ESTIMATED STATE
    %estimated_deputy_rel_LVLH = AbsECI2RelHill_V1_1(estimated_satellites_iner_ECI(:,i+1), mi);
    %[u_estimated, temp, estimated_params] = ControlEvaluation_V1_4(estimated_deputy_rel_LVLH, estimated_satellites_iner_ECI, t, params);
    %UestimatedStory(:,i) = [u_estimated(1,1:end)';u_estimated(2,1:end)';u_estimated(3,1:end)'];
    %params.u = memo_params_u;
    
    ObserverTest.satellites_iner_ECI_estimated = ObserverTest.estimated_satellites_iner_ECI(:,i+1);
    ObserverTest.deputy_rel_LVLH_estimated = AbsECI2RelHill_V1_2(ObserverTest.satellites_iner_ECI_estimated, params.mi);
    satellites_iner_ECI_estimated = ObserverTest.satellites_iner_ECI_estimated;
    deputy_rel_LVLH_estimated = ObserverTest.deputy_rel_LVLH_estimated; 
    
end
%% END POSITION OBSERVER


%%ATTITUDE OBSERVER

if(ObserverTest.AttitudeObserverOn==1)
    
    for k=1:ObserverTest.Nagents,
        
        tspan_att = tspan; %THEY SHOULD BE THE SAME
        
        %Attitude_Observer_V2_1(satellites_attitude, params,i,tspan_att,k);
        Attitude_ObserverOpt_V2_1(satellites_attitude, params,i,tspan_att,k);
        
        ObserverTest.estimated_satellites_attitude(1+7*(k-1):7+7*(k-1),i+1) = Agent(k).attitude_xHatUKF(1:7,i+1); 
        
    end
    ObserverTest.satellites_attitude_estimated = ObserverTest.estimated_satellites_attitude(:,i+1);
    satellites_attitude_estimated = ObserverTest.satellites_attitude_estimated;
end


%position_error = satellites_iner_ECI - satellites_iner_ECI_estimated
%attitude_error = satellites_attitude - satellites_attitude_estimated





        