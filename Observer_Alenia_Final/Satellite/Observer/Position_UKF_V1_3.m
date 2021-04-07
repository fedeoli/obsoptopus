%UKF to estimate position and speeds of the agents

for k=1:ObserverTest.Nagents,
    
    
    if(Agent(k).UWBcorrectionActive(i+1)==1)
        myGPS = Agent(k).GPSopt(:,i+1);
    else
        myGPS = Agent(k).GPS(:,i+1);
    end
    
    %RESET OF P after a given amount of samples
    if(mod(i+1,ObserverTest.position_P_reset_aftersamples)==0)
        %disp(['P Reset - det(P)=  '   num2str(det(Agent(k).P)) ' substituted with det(Preset) = ' num2str(det(reshape(ObserverTest.PiAll(Ntest,:),size(ObserverTest.Pi))))])
        Agent(k).P = reshape(ObserverTest.PiAll(ObserverTest.Ntest,:),size(ObserverTest.Pi));
    end
    
    ObserverTest.Pa = blkdiag(Agent(k).P,Agent(k).Q);
    [ObserverTest.W,ObserverTest.CSI,Psat] = sigmaPoints_sat(Agent(k).xaUKF,ObserverTest.Pa,params);
    ObserverTest.Pa = Psat;
    %Sigma point evolution
    ObserverTest.CSI_X_meno = 0*ObserverTest.CSI(1:6,:);
    for j=1:2*ObserverTest.Na+1,
        ObserverTest.CSI_jth = ObserverTest.CSI(:,j)';
        %DX =  InertialDynamicsIntegratorV11_EstDecentr1(CSI(1:6,j), params,k,CSI(6+1:end,j),control_time_index);
        %ObserverTest.CSI_X_meno(:,j) = CSI(1:6,j) + DX*time_step ;
        X = rk4_V1_1_decentralized(@InertialDynamicsIntegrator_V2_1_decentralized, tspan,ObserverTest.CSI_jth(1:6), params,k);
        %X(:,end) = X(:,end) + [ObserverTest.CSI(7:9,j);ObserverTest.CSI(7:9,j)]*((time_step)^2/2 + (time_step));
        ObserverTest.CSI_X_meno(:,j) = X(:,end);
    end
    
    % stima valor medio
    xHatMeno = ObserverTest.CSI_X_meno*ObserverTest.W;
    Agent(k).xHatMeno = xHatMeno;
    % stima a priori
    
    % Aggiornamento matrice di covarianza
    MatCOV = Agent(k).P;
    for j=1:2*ObserverTest.Na+1,
        MatCOV = MatCOV + ObserverTest.W(j)*(ObserverTest.CSI_X_meno(:,j)-xHatMeno)*(ObserverTest.CSI_X_meno(:,j)-xHatMeno)';
    end
    Pmeno = MatCOV;
    Agent(k).Pmeno = Pmeno;
    
    % GPS + UWB attesi per ogni sigma point
    % ZETA: the first three are the three coordinates of the GPS while
    % the other Nagents-1 are the distance measured with the other
    % ones (the own distance is cleared off).
    for j=1:2*ObserverTest.Na+1,
        if(ObserverTest.AddDerivativeEstimation ~= 0)% && i >ObserverTest.PseudoDerivative_d)
            ZETA(:,j) = [ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS,j);ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS+3,j)];
        else
            ZETA(:,j) = [ObserverTest.CSI_X_meno(ObserverTest.ExtractGPS,j)];
        end
    end
    % expected value of the measures (GPS) and speed
    zHatMeno = ZETA*ObserverTest.W;
    
    % matrice covarianza Pzeta
    MatCOV = zeros(length(ZETA(:,1)),length(ZETA(:,1)));
    for j=1:2*ObserverTest.Na+1,
        MatCOV = MatCOV + ObserverTest.W(j)*(ZETA(:,j)-zHatMeno)*(ZETA(:,j)-zHatMeno)';
    end
    Pzeta = MatCOV + Agent(k).R;
    
    % matrice covarianza mista Pxz
    MatCOV = zeros(length(ObserverTest.CSI_X_meno(:,1)),length(ZETA(:,1)));
    for j=1:2*ObserverTest.Na+1,
        MatCOV = MatCOV + ObserverTest.W(j)*(ObserverTest.CSI_X_meno(:,j)-xHatMeno)*(ZETA(:,j)-zHatMeno)';
    end
    Pxz = MatCOV;
    
    % innovation ?
    if(ObserverTest.AddDerivativeEstimation ~= 0)% && i >ObserverTest.PseudoDerivative_d))
        if(ObserverTest.UWB_StepsBackInTime > 0 && ObserverTest.EnableSpeedOptimizationWhenSteppingBack)%then we also optimize the speed and we use it
            innovation = [myGPS; Agent(k).GPSpeedOpt(:,i+1)] - zHatMeno;
        else
            innovation = [myGPS; Agent(k).GPSderivative(:,i+1)] - zHatMeno;
        end
    else
        innovation = [myGPS] - zHatMeno;
    end
    
    % guadagno UKF
    K = ObserverTest.CorrectionON*Pxz*pinv(Pzeta);
    
    %This is the corrected estimation at time i+1-Sensor_delay
    xHat = Agent(k).xHatMeno + K*innovation;

    
    if(iteration < ObserverTest.OptimizationsEachSampling ) %it is not the last that need to update everything
        Agent(k).xHatUKFmultiple = xHat;
    else
        % Aggiornamento matrice di covarianza at time i+1-Sensor_delay
        P = Agent(k).Pmeno - K*Pzeta*K';
        
        %prediction of the estimate
        Agent(k).xHat = xHat; %the first component is the variations with respect to r
        Agent(k).xHatUKF(:,i+1) = Agent(k).xHat;
        Agent(k).xaUKF = [Agent(k).xHat; zeros(ObserverTest.Ndisturbance,1)]; %[xHat; zeros(Ndisturbance,1)]; %note that xaUKF contains the estimated state at time i+1-Sensor_delay
        Agent(k).P = P;
        
        Agent(k).iner_ECI_EstimationError(:,i+1) = Agent(k).iner_ECI(:,i+1) - Agent(k).xHatUKF(:,i+1);
        Agent(k).NormEstimationErrorXYZ(i+1) = norm(Agent(k).iner_ECI_EstimationError(1:3,i+1));
        %XhatEuler(1+(k-1)*6:6+(k-1)*6,i+1) = Agent(k).xHatUKF(:,i+1);
        ObserverTest.estimated_satellites_iner_ECI(1+(k-1)*6:6+(k-1)*6,i+1) = Agent(k).xHatUKF(:,i+1);
    end
end