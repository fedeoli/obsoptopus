
SkipChiefDidAll = 0;
k = 1;
while( k <= ObserverTest.Nagents && SkipChiefDidAll ==0);
    
    
    %% GPS OPTIMIZATION VIA UWB
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
    
    myGPS_obs = myGPS;
    ObserverTest.Kth_agent = k;
    ObserverTest.actual_time_index = i+1;
    ObserverTest.u_time_index = i;
    
    %Some of the distances might not be received
    Agent(k).SuccessfullyReceivedData(1:end-1,:) = Agent(k).SuccessfullyReceivedData(2:end,:);
    Agent(k).SuccessfullyReceivedData(end,:) = ones(1,ObserverTest.Nagents);
    if(ObserverTest.UWBDropMessages == 1)
        Agent(k).SuccessfullyReceivedData(end,:) = SendMessagesMatrix_v1(Agent(k).SuccessfullyReceivedData(end,:),ObserverTest.UWBDropMessagesP,ObserverTest.UWBDropMessagesR);
    end
    
    %check if GPS optimization has to be done
    if((ObserverTest.StartSoonUWB == 1 || ((norm(myGPS_obs - Agent(k).xHatUKF(1:3,i+1)) < ObserverTest.UWBonGPSerror && i>ObserverTest.UWBOptimizationNoBeforeThan)) && ObserverTest.UWBoptimizationOn) )
        %do not start yet otimization from past estimated state
        
        %SETTING THE INITIAL CONDITIONS FOR OPTIMIZATION
        %use GPS
        if( ObserverTest.UWBinitialConditionOptimizationIsGPS==1)
            ObserverTest.X0_start_Kth_agent = myGPS;% Agent(k).GPSoffsetEstimated(:,i)];
            if(ObserverTest.CentralizedOptimization == 1)
                for jj=1:Nagents,
                    ObserverTest.X0_start_Kth_agent([1:3] + 3*(jj-1)) = Agent(jj).GPS(:,i+1);
                end
            end
            
        else
            %use the previous estimations
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
        
        Agent(k).UWBcorrectionActive(i+1) = 1;
        
        %%%%%%%%%%%%% GPS OPTIMIZATION %%%%%%%%%%%%% 
        % Get apriori estimate 
        Chi = reshape(ObserverTest.APrioriEstimationXYZ,3,4)';
        
        % Get relative distances
        adjmat_UWB = ObserverTest.AprioriAdjacencyMatrix;
        
        % Get GPS measurements
        GPS = zeros(ObserverTest.Nagents,3);
        for z = 1:ObserverTest.Nagents
            GPS(z,:) = Agent(z).GPS(:,i);
        end
        
        % weights selection - geometric method
        %%% SIGMA OPTIMIZED %%%
        theta = 0.05;
        beta = 0.3;
        
        % optimize GPS
        opt = Position_opt_cloud_num_v5(Chi, GPS, adjmat_UWB, theta, beta);
        NewGPS = reshape(opt.Chi_est,1,3*ObserverTest.Nagents);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
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
        
        %Agent(k).GPSoffsetEstimated(:,i+1) = NewGPS(4:6);
    end
    k = k + 1;
end
