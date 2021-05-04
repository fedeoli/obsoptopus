 function  IterativePseudoDerivativeCodeAgentGPS(k_th,k_agent)   
   %Calcolo del rapporto incrementale mediando i primi e gli utlimi c campioni a distanza d 
global Agent ObserverTest

% if(k_th>=2*Test.PseudoDerivative_d)
%     disp('ciao');
% end

for j=1:3,
    
    j1 = max(1,k_th-ObserverTest.PseudoDerivative_c+1);
    j2 = max(1,k_th-ObserverTest.PseudoDerivative_d);
    j3 = max(1,k_th-ObserverTest.PseudoDerivative_d+ObserverTest.PseudoDerivative_c-1);
    if(k_th > ObserverTest.FilterSteps2Regime2EvaluateDerivatives && ObserverTest.GPSfilterOn==1)
        temp1 = mean(Agent(k_agent).GPSfiltered(j,j1:k_th));
        temp2 = mean(Agent(k_agent).GPSfiltered(j,j2:j3));
    else
        temp1 = mean(Agent(k_agent).GPS(j,j1:k_th));
        temp2 = mean(Agent(k_agent).GPS(j,j2:j3));
    end
    
    if(ObserverTest.ZeroErrors == 1) %no error
        Agent(k_agent).GPSderivative(j, k_th) = Agent(k_agent).iner_ECI( 3+j,k_th);
    else
        Agent(k_agent).GPSderivative(j, k_th) = (temp1-temp2) / (ObserverTest.time_step*max(1,k_th-j2-(j3-j2)));  
    end
    %temp1 = mean(Agent(k_agent).iner_ECI(j,j1:k_th));
    %temp2 = mean(Agent(k_agent).iner_ECI(j,j2:j3));
    %Agent(k_agent).GPSderivative(j, k_th) = (Agent(k_agent).iner_ECI( j,k_th)- Agent(k_agent).iner_ECI( j,max(1,k_th-1))) /ObserverTest.SamplingTime;  
    
end



