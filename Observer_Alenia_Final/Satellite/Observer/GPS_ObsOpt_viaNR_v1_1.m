function x = GPS_ObsOpt_viaNR_v1_1(x)

global ObserverTest Agent params


WGPS = ObserverTest.Weight_GPS_back(1),
WSIGMA = 
WUWB







    
    %% DECENTRALIZED
    x_position = x(1:3);
    for n=0:ObserverTest.UWB_StepsBackInTime,
        if(n==0)
            chi_i = x_position;
            if(ObserverTest.UWB_StepsBackInTime > 0)
                if(ObserverTest.EnableSpeedOptimizationWhenSteppingBack == 1)
                    temp_chi_i = x(1:6);
                else
                    temp_chi_i = [x_position; Agent(ObserverTest.Kth_agent).xHatUKF(4:6,ObserverTest.actual_time_index) ];
                end
            end
        else %we simulate back in time
            %chi_i = x_position - ObserverTest.time_step*sum(Agent(ObserverTest.Kth_agent).GPSpeed(:,ObserverTest.actual_time_index-n+1:ObserverTest.actual_time_index),2);
            %Time range of integration bewteen samples
            tspan = ObserverTest.time_step*(ObserverTest.actual_time_index - (n)*ObserverTest.WindowStepMultiplier): ObserverTest.time_step*ObserverTest.WindowStepMultiplier: ObserverTest.time_step*(ObserverTest.actual_time_index - (n-1)*ObserverTest.WindowStepMultiplier);
            X0 = temp_chi_i;
            for jj =1:ObserverTest.WindowStepMultiplier,
                ObserverTest.u_time_index = ObserverTest.actual_time_index - (n-1)*ObserverTest.WindowStepMultiplier - jj; %the u index is the end -1 of the time range that is considered
                %make reverse jump before integrating as it is done after
                %integration when forward in time
                X = rk4_V1_1_decentralized(@InertialDynamicsIntegrator_V1_1_decentralized_back, tspan(end-jj:end-jj+1), X0, params,ObserverTest.Kth_agent);
                X0 = X(:,end);
            end
            chi_i = X(1:3,end);
            temp_chi_i = X(1:6,end);
        end
        for j=1:ObserverTest.Nagents, %Evaluate the squared distance among agent Kth_agent and the others, if received
            if(Agent(ObserverTest.Kth_agent).SuccessfullyReceivedData(end-n*ObserverTest.WindowStepMultiplier,j)==1  && ObserverTest.Kth_agent~=j)%if the data of j-th agent is arrived
                chi_j = Agent(ObserverTest.Kth_agent).GotData(ObserverTest.WindowStepMultiplier*(ObserverTest.UWB_StepsBackInTime+1-n), [1:3]+3*(j-1))';
                temp =  abs(norm(chi_i - chi_j) -  Agent(ObserverTest.Kth_agent).MeasuredDistances(ObserverTest.actual_time_index-ObserverTest.WindowStepMultiplier*n,j)) ;
                if( temp < ObserverTest.DeadZoneUWB)
                    temp = 0;
                end
                %set to zero if data from j-th agent has been lost
                J1 = J1 +( ObserverTest.Weight_UWB_back(n+1)/((ObserverTest.Nagents-1)*(ObserverTest.UWB_StepsBackInTime+1)))*(temp*ObserverTest.MeasureMagnifier).^ObserverTest.expJUWB;
            end
        end
        
        if(Agent(ObserverTest.Kth_agent).SuccessfullyReadGPS(end-n*ObserverTest.WindowStepMultiplier))

%             %OBSOLETE            
%             if(ObserverTest.GPSoffset)
%                 GPSoffsetEstimated = x(4:6);
%             else
%                 GPSoffsetEstimated = 0*x_position;
%             end
            GPSoffsetEstimated = 0*x_position;
            
            GPS = Agent(ObserverTest.Kth_agent).GPS(:,ObserverTest.actual_time_index-n*ObserverTest.WindowStepMultiplier);
            if(ObserverTest.GPSfilterOn == 1)
                GPS = Agent(ObserverTest.Kth_agent).GPSfiltered(:,ObserverTest.actual_time_index-n*ObserverTest.WindowStepMultiplier);
            end
            temp = norm(GPS + GPSoffsetEstimated - chi_i);
            if( temp < ObserverTest.DeadZoneGPS)
                temp = 0;
            end
            
            J2 = J2 + (ObserverTest.Weight_GPS_back(n+1)/(ObserverTest.UWB_StepsBackInTime+1))*(temp*ObserverTest.MeasureMagnifier)^ObserverTest.expJGPS;
            %accounting for speed error
            if(ObserverTest.EnableSpeedOptimizationWhenSteppingBack == 1)
                J4 = J4 + (ObserverTest.Weight_GPSpeed/(ObserverTest.UWB_StepsBackInTime+1))*(ObserverTest.MeasureMagnifier*norm(temp_chi_i(4:6)-Agent(ObserverTest.Kth_agent).GPSpeed(:,ObserverTest.actual_time_index-n*ObserverTest.WindowStepMultiplier)))^ObserverTest.expJGPSpeed;
            end
        end
    end
    
    %faster if the estimate is far from GPS error bound (we drop the phase lag for faster convergence)
    WeightSigma = ObserverTest.Weight_Innovation;
    if(Agent(ObserverTest.Kth_agent).MismatchAprioriDistances(ObserverTest.actual_time_index) > ObserverTest.AprioriMismatch4FastConvergence)
       WeightSigma = ObserverTest.Weight_Innovation_Fast;
    end
    if(mod(ObserverTest.Weight_Innovation_Reset_aftersamples,ObserverTest.actual_time_index)==0)
        WeightSigma = ObserverTest.Weight_Innovation_Fast_Reset;
    end
    J3 = WeightSigma*(norm(x_position(1:3) - ObserverTest.X0_start_Kth_agent(1:3))*ObserverTest.MeasureMagnifier)^ObserverTest.expJSIGMA;
    
    if(ObserverTest.EnableSpeedOptimizationWhenSteppingBack == 1)
        J5 = ObserverTest.Weight_GPSpeedSigma*(ObserverTest.MeasureMagnifier*norm( x(4:6) - ObserverTest.X0_start_Kth_agent(4:6) )  )^ObserverTest.expJGPSpeedSigma;
    end
end

J =  J1 + J2 + J3 + J4 + J5;

if(ObserverTest.OptimizationCycles.count==0)
    ObserverTest.OptimizationCycles.Initial = J;
end
ObserverTest.OptimizationCycles.count = ObserverTest.OptimizationCycles.count +1;
ObserverTest.OptimizationCycles.End = J;


%J2 = Weight_GPS*norm(myGPS_obs - GPSoffsetEstimated - x_position)^4;
%ObserverTest.N_WsigmaOffset*(norm(GPSoffsetEstimated))^4;%-X0_start_Kth_agent(4:6))^4);
%J = Weight_UWB*sum(abs(AdmRowK_opt' - Hat_AdmRowK_opt )) + ...
%                         Weight_GPS*norm(myGPS_obs - x)^2 + ...
%                        Weight_Innovation_Fast*norm(x - X0_start_Kth_agent)^2;
