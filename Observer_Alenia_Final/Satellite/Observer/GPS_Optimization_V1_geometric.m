%% GPS OTPIMIZATON - GEOMETRIC METHOD

global ObserverTest Agent

% flag and iteration init 
SkipChiefDidAll = 0;
k = 1;
ObserverTest.UWBoptimizationOn = 1;

% loop over the agents
while( k <= ObserverTest.Nagents && SkipChiefDidAll ==0);
        
%     GPS OPTIMIZATION VIA UWB

    % get current GPS
    if(iteration == 1)
        myGPS = Agent(k).GPS(:,i+1);
        if(ObserverTest.GPSfilterOn == 1)
            myGPS = Agent(k).GPSfiltered(:,i+1);
        end
    else
        myGPS = Agent(k).GPSopt(:,i+1); %previous optimized GPS is considered
        for jj=1:ObserverTest.Nagents,
            if(ObserverTest.TrasmittedTruePositions == 1)
                Agent(n).GotData(end,(jj-1)*3+1:(jj-1)*3+3) =  Agent(jj).iner_ECI(1:3,i+1);
            else
                Agent(n).GotData(end,(jj-1)*3+1:(jj-1)*3+3) =  Agent(jj).xHatUKFmultiple(1:3);%NOTE HERE, not the standard xHatUKF is shared but the xhat evaluate by multiple iterations on the same sample
            end
        end
    end
    
    % iteration data init
    myGPS_obs = myGPS;
    ObserverTest.Kth_agent = k;
    ObserverTest.actual_time_index = i+1;
    ObserverTest.u_time_index = i;
    
    % Some of the distances might not be received
    Agent(k).SuccessfullyReceivedData(1:end-1,:) = Agent(k).SuccessfullyReceivedData(2:end,:);
    Agent(k).SuccessfullyReceivedData(end,:) = ones(1,ObserverTest.Nagents);
    if(ObserverTest.UWBDropMessages == 1)
        Agent(k).SuccessfullyReceivedData(end,:) = SendMessagesMatrix_v1(Agent(k).SuccessfullyReceivedData(end,:),ObserverTest.UWBDropMessagesP,ObserverTest.UWBDropMessagesR);
    end
    
    % flag for GPS optimization
    GPS_flag = ((ObserverTest.StartSoonUWB == 1 || ((norm(myGPS_obs - Agent(k).xHatUKF(1:3,i+1)) < ObserverTest.UWBonGPSerror && i>ObserverTest.UWBOptimizationNoBeforeThan)) && ObserverTest.UWBoptimizationOn) );

    % check if GPS optimization has to be done
    if( GPS_flag )
        % do not start yet otimization from past estimated state
        
        % SETTING THE INITIAL CONDITIONS FOR OPTIMIZATION
        % use GPS
        if( ObserverTest.UWBinitialConditionOptimizationIsGPS==1)
            ObserverTest.X0_start_Kth_agent = myGPS;% Agent(k).GPSoffsetEstimated(:,i)];
            if(ObserverTest.CentralizedOptimization == 1)
                for jj=1:Nagents,
                    ObserverTest.X0_start_Kth_agent([1:3] + 3*(jj-1)) = Agent(jj).GPS(:,i+1);
                end
            end
            
        else
            % use the previous estimations
            if(iteration > 1)
                if(ObserverTest.CentralizedOptimization == 1)
                    for jj=1:Nagents,
                        ObserverTest.X0_start_Kth_agent([1:3] + 3*(jj-1)) = Agent(k).xHatUKFmultiple(1:3)';
                    end
                else
                    ObserverTest.X0_start_Kth_agent = Agent(k).xHatUKFmultiple(1:3)';
                end
                
            else
                if(ObserverTest.CentralizedOptimization == 1)
                    for jj=1:ObserverTest.Nagents,
                        ObserverTest.X0_start_Kth_agent(1 + 3*(jj-1),1) = Agent(jj).xHatUKF(1,i+1);
                        ObserverTest.X0_start_Kth_agent(2 + 3*(jj-1),1) = Agent(jj).xHatUKF(2,i+1);
                        ObserverTest.X0_start_Kth_agent(3 + 3*(jj-1),1) = Agent(jj).xHatUKF(3,i+1);
                        if(ObserverTest.UWB_StepsBackInTime > 0) %in this case also the speed is estimated propagating back in time
                            ObserverTest.X0_start_Kth_agent(4+ 3*(jj-1),1) = Agent(jj).xHatUKF(4,i+1);
                            ObserverTest.X0_start_Kth_agent(5+ 3*(jj-1),1) = Agent(jj).xHatUKF(5,i+1);
                            ObserverTest.X0_start_Kth_agent(6+ 3*(jj-1),1) = Agent(jj).xHatUKF(6,i+1);
                        end
                    end
                else
                    ObserverTest.X0_start_Kth_agent(1,1) = Agent(k).xHatUKF(1,i+1);
                    ObserverTest.X0_start_Kth_agent(2,1) = Agent(k).xHatUKF(2,i+1);
                    ObserverTest.X0_start_Kth_agent(3,1) = Agent(k).xHatUKF(3,i+1);
                    if(ObserverTest.UWB_StepsBackInTime > 0) %in this case also the speed is estimated propagating back in time
                        ObserverTest.X0_start_Kth_agent(4,1) = Agent(k).xHatUKF(4,i+1);
                        ObserverTest.X0_start_Kth_agent(5,1) = Agent(k).xHatUKF(5,i+1);
                        ObserverTest.X0_start_Kth_agent(6,1) = Agent(k).xHatUKF(6,i+1);
                    end
                        
                end
                
            end
            
        end
        
        % set GPS optimization flag for the current iteration
        Agent(k).UWBcorrectionActive(i+1) = 1;
        
        %%%%%%%%%%%%% GPS OPTIMIZATION %%%%%%%%%%%%% 
        % Get apriori estimate 
        Chi = reshape(ObserverTest.APrioriEstimationXYZ,3,ObserverTest.Nagents)';
        
        % Get past estimation - all agents
        Chi_past = zeros(ObserverTest.Nagents,3);
        for z = 1:ObserverTest.Nagents
            Chi_past(z,:) = Agent(z).xHatUKF(1:3,i);
        end
        
        % Get relative distances
        adjmat_UWB = ObserverTest.MeasuredDistances;
        
        % Get GPS measurements
        if strcmp(ObserverTest.projection,'Chi')
            GPS = reshape(Agent(k).GPS(:,i+1),1,3);
        elseif strcmp(ObserverTest.projection,'GPS')
            GPS = zeros(ObserverTest.Nagents,3);
            for z = 1:ObserverTest.Nagents
                GPS(z,:) = reshape(Agent(z).GPS(:,i+1),1,3);
            end
        end
        
        % optimize GPS
%         opt = Position_opt_cloud_num_v6_dec(Chi, Chi_past, GPS, adjmat_UWB, ObserverTest.Kth_agent, ObserverTest.theta, ObserverTest.beta, ObserverTest.check_distance);
        opt = Position_opt_cloud_num_v5_dec(Chi, GPS, adjmat_UWB, ObserverTest.Kth_agent, ObserverTest.theta, ObserverTest.beta, ObserverTest.check_distance,...
                ObserverTest.projection);
        NewGPS = [reshape(opt.Chi_est,1,3), zeros(1,3)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % post optimization assignments
        if(ObserverTest.CentralizedOptimization==1)
            myGPS = NewGPS([1:3]+(k-1)*3);
            Agent(k).GPSopt(:,i+1) = NewGPS([1:3]+(k-1)*3);
            %in case of CentralizedOptimization, only one
            %perform it, the Chief, and then pass it to
            %all the others the optimized GPS for all
            for kk=1:ObserverTest.Nagents,
                Agent(kk).GPSopt(:,i+1) = NewGPS([1:3]+(kk-1)*3);
                Agent(kk).UWBcorrectionActive(i+1) = 1;
                if(ObserverTest.UWB_StepsBackInTime > 0 && ObserverTest.EnableSpeedOptimizationWhenSteppingBack)
                    Agent(kk).GPSpeedOpt(:,i+1) = NewGPS([4:6]+(kk-1)*3);
                end
            end
            SkipChiefDidAll = 1;%NOTE THAT K IS CHANGED WITHIN THE FOR TO BREAK IT
        else
            myGPS = NewGPS(1:3);
            Agent(k).GPSopt(:,i+1) = myGPS;
            if(ObserverTest.UWB_StepsBackInTime > 0 && ObserverTest.EnableSpeedOptimizationWhenSteppingBack)
                Agent(k).GPSpeedOpt(:,i+1) = NewGPS(4:6);
            end
        end
        
        % compute estimation error
        Agent(k).iner_ECI_EstimationError(:,i+1) = Agent(k).iner_ECI(:,i+1) - [Agent(k).GPSopt(:,i+1); Agent(k).GPSpeedOpt(:,i+1)];
        Agent(k).NormEstimationErrorXYZ(i+1) = norm(Agent(k).iner_ECI_EstimationError(1:3,i+1));
        
        % set new GPS
        ObserverTest.estimated_satellites_iner_ECI(1+(k-1)*6:3+(k-1)*6,i+1) = Agent(k).GPSopt(:,i+1);
        ObserverTest.estimated_satellites_iner_ECI(4+(k-1)*6:6+(k-1)*6,i+1) = Agent(k).GPSpeedOpt(:,i+1);
        
    else  
        Agent(k).GPSopt(:,i+1) = Agent(k).xHatUKF(1:3,i+1);
        Agent(k).GPSpeedOpt(:,i+1) = Agent(k).xHatUKF(4:6,i+1);

        % set new GPS        
        ObserverTest.estimated_satellites_iner_ECI(1+(k-1)*6:3+(k-1)*6,i+1) = Agent(k).GPSopt(:,i+1);
        ObserverTest.estimated_satellites_iner_ECI(4+(k-1)*6:6+(k-1)*6,i+1) = Agent(k).GPSpeedOpt(:,i+1);
        
        % compute estimation error
        Agent(k).iner_ECI_EstimationError(:,i+1) = Agent(k).iner_ECI(:,i+1) - [Agent(k).GPSopt(:,i+1); Agent(k).GPSpeedOpt(:,i+1)];
        Agent(k).NormEstimationErrorXYZ(i+1) = norm(Agent(k).iner_ECI_EstimationError(1:3,i+1));
        
        
    end
    
    % update agent ID
    k = k + 1;
end

