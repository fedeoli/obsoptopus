function J = attitudeObsOptimization_V1_1(x)

global ObserverTest Agent params

J1 = 0;
J2 = 0;
J3 = 0;
    %% DECENTRALIZED
%ObserverTest.attitudeBackSamplesWindow = 10; %how many samples to fit (window size for data matching)
%ObserverTest.attitudeBackIntersample = 1;%number of steps among samples

Window = ObserverTest.attitudeBackSamplesWindow;
Intersamples = ObserverTest.attitudeBackIntersample;
TimeIndex = ObserverTest.attitudeTimeIndex;
TimeStep = ObserverTest.time_step;

state = x; %state at present time ObserverTest.attitudeTimeIndex
for n=0:((Window-1)*(Intersamples)),
    
    %chi_i = x_position - ObserverTest.time_step*sum(Agent(ObserverTest.Kth_agent).GPSpeed(:,ObserverTest.actual_time_index-n+1:ObserverTest.actual_time_index),2);
     %Time range of integration bewteen samples
     if(n>0)
        tspan = [(TimeIndex-n):1:(TimeIndex-n+1)]*TimeStep;
        %make reverse jump before integrating since it is done after integration when forward in time
        ObserverTest.attitude_u_backIntime_index = (TimeIndex-n);
        
        X = rk4_V1_1_decentralized(@AttitudeDynamicsOpt_V1_1_decentralized_back, tspan, state, params,ObserverTest.Kth_agent);
        state = X(:,end);
     end
     
     %we evaluate the mismatch between the estimated measures and the
     %sampled ones each ObserverTest.attitudeBackIntersample.
     if(0 == mod(n,(ObserverTest.attitudeBackIntersample-1)))
        %evaluate the measure associated with "state", the optimization variable
        [MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(Agent(ObserverTest.Kth_agent).iner_ECI(1:3,ObserverTest.attitudeTimeIndex-n)',...
        state(1:4)', state(5:7)',state(8:10)',state(11:13)',1); 
        J1 = J1 + norm(Agent(ObserverTest.Kth_agent).magneto(:, ObserverTest.attitudeTimeIndex-n) - MagTemp)^2; %weight on magnetic mismatch
        J2 = J2 + norm(Agent(ObserverTest.Kth_agent).gyros(:, ObserverTest.attitudeTimeIndex-n) - GyrosTemp)^2; %weight on omega mismatch
        J3 = ObserverTest.attitudeWnormquaternion*abs(quatnorm(state(1:4)')-1);%weight on quaternion constraint
     end  
    
end

J =  J1 + J2 + J3 ;

disp(['J: ' num2str(J)]);

ObserverTest.attitudeOptimizationCounter = ObserverTest.attitudeOptimizationCounter + 1;
if(ObserverTest.attitudeOptimizationCounter==1)
    ObserverTest.attitudeOptimizationCyclesInitial = J;
else
    ObserverTest.attitudeOptimizationCyclesEnd = J;
end
    

