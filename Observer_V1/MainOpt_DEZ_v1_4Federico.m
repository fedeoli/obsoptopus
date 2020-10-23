% Main Optimization
clc
clear all
close all

global DynOpt

DynOpt.Ts = 0.0005; %: Sampling Time
simulationModel = 1;


%% PLANT data
Ki = 1.5; %good 1.5*2;
Kp = 0.22; %good 0.22*30;
Kd = 0.0; %good 1.0;
beta = 140;%good 30;
u0 = 1.0;

%model simulation
DynOpt.c0 = 45;
DynOpt.c1 = 20;
DynOpt.gamma0 = 10;
DynOpt.gamma1 = 0.01;
DynOpt.beta = beta;
DynOpt.StateDim = 4;
DynOpt.Zh = 1.5;
DynOpt.Zf = 1.3;
DynOpt.Ki = Ki*1;
DynOpt.Kp = Kp*1;
DynOpt.Kd = Kd*1;
DynOpt.beta = beta;
DynOpt.u0 = u0;
DynOpt.gainPHD = 1/(1.5e-4);

A = [0, 1,0,0; 0,0,1,0; ...
        0,-DynOpt.c0+DynOpt.gamma0,-DynOpt.c1,DynOpt.gamma1;
        -DynOpt.Ki*DynOpt.gainPHD*DynOpt.beta,-DynOpt.Kp*DynOpt.gainPHD*DynOpt.beta,-DynOpt.Kd*DynOpt.gainPHD*DynOpt.beta,-DynOpt.beta];
 eig(A)
 eig(eye(size(A))+DynOpt.Ts*A)
    
    if(simulationModel == 1)
        time = [0:DynOpt.Ts:0.3];
        N = length(time);
        stateStory = zeros(4,N);
        stateStory(:,1) = [0;0.01;0;0];
        DynOpt.U = zeros(2,N);

        
        stop = 0;
        for k=2:N,
            if(stop==0)
                DynOpt.U(2,k-1) = 1+0.2*sin(2*pi*1*time(k-1));
                DynOpt.U(1,k-1) = -DynOpt.Ki*DynOpt.gainPHD*stateStory(1,k-1)-DynOpt.Kp*DynOpt.gainPHD*stateStory(2,k-1)-DynOpt.Kd*DynOpt.gainPHD*stateStory(3,k-1);
                stateStory(:,k) = PlantJumpMap(stateStory(:,k-1),k-1,1,0);
                if(stateStory(2,k)>=0.5)
                    stop = 1;
                end
            end
        end
        
        figure(1)
        ax(1)=subplot(4,1,1);
        plot(time,stateStory(1,:),'LineWidth',2);
        grid on
        ylabel('$\int z$')
        title('Simulation test ')
        ax(2)=subplot(4,1,2);
        plot(time,stateStory(2,:),'LineWidth',2);
        grid on
        ylabel('z')
        ax(3)=subplot(4,1,3);
        plot(time,stateStory(3,:),'LineWidth',2);
        grid on
        ylabel('$\dot z$')
        ax(4)=subplot(4,1,4);
        plot(time,stateStory(4,:),'-',time,max(0,(DynOpt.U(2,:)-1))*max(stateStory(4,:)),'--',time,DynOpt.U(1,:),':','LineWidth',2);
        grid on
        ylabel('u_m')
        linkaxes(ax,'x');
        
        
    end
%%


if(simulationModel)
    z = stateStory(2,:);
    u_c = DynOpt.U(1,:);
    u_m = DynOpt.U(1,:);
    elong = DynOpt.U(2,:);
    state = stateStory;
else
    state = [int_z;z;dz;u_m_est];
end

DynOpt.c1_derivative = 6;
DynOpt.d1_derivative = 20;%has to be smaller than DynOpt.Nts
dz = (mypseudo_derivative(time,z,DynOpt.c1_derivative,DynOpt.d1_derivative))';
DynOpt.c1_dderivative = 8;
DynOpt.d1_dderivative = 30; %has to be smaller than DynOpt.Nts
ddz = (mypseudo_derivative(time,dz,DynOpt.c1_dderivative,DynOpt.d1_dderivative))';
int_z = DynOpt.Ts^(-1)*cumtrapz(time,z);
u_m_est = u0*lsim(c2d(tf(1,[1/beta 1]),DynOpt.Ts,'zoh'),u_c,time)';





%% OBSERVER

DynOpt.Nts = 50 % Inter-sampling Time, number of samples within intersampled data
DynOpt.w = 8 %number of inter-sampled data in a window: the total number of sampled data in a window are (DynOpt.w-1)*DynOpt.Nts+1
DynOpt.Acon = eye(6); DynOpt.Acon(5,5) = -1*0;DynOpt.Acon(6,6) = -1*0;
DynOpt.Bcon = zeros(6,1);
DynOpt.Aeq = zeros(6);
DynOpt.Beq = zeros(6,1);
DynOpt.lb = [-5E1,-1E2,-4E2,-5E1,2,0];
DynOpt.ub = [5E1,1E2,4E2,5E1,30,30];
DynOpt.Weight  = ones(1,DynOpt.w); %window weight size, last on the more recent t_k
DynOpt.dWeight  = 0.0*ones(1,DynOpt.w); %window weight size, last on the more recent t_k
DynOpt.dWeight(1) = 0; DynOpt.dWeight(end) = 0;
DynOpt.ddWeight  = 0*ones(1,DynOpt.w); %window weight size, last on the more recent t_k
DynOpt.ddWeight(1) = 0; DynOpt.ddWeight(end) = 0;
%DynOpt.Xtrue = [int_z(1);z(1);dz(1);u_m(1);DynOpt.gamma0;DynOpt.gamma1]; %true state and parameters
DynOpt.Xtrue = [state(:,1);DynOpt.gamma0;DynOpt.gamma1]; 
DynOpt.X  = 1.2*DynOpt.Xtrue;
DynOpt.WindowSamples = max(2,DynOpt.Nts*(DynOpt.w-1)+1);
DynOpt.U =  zeros(2,DynOpt.WindowSamples); %+1 since we need to maintain the previous time for integration.
DynOpt.Y =  zeros(3,DynOpt.w);
DynOpt.StateDim = 4;
myoptioptions = optimoptions(@fminunc,'Algorithm','quasi-newton',  'MaxIter', 50,'display','off'); % 'Algorithm', 'trust-region','SpecifyObjectiveGradient', true,
%myoptioptions = optimoptions(@fmincon,  'MaxIter', 100,'display','off'); % 'Algorithm', 'trust-region','SpecifyObjectiveGradient', true,


DynOpt.TestDynamics = 0; %different from zero to test only the backward integration and see if it is correct: please set the correct value in the optimized variable initialization
DynOpt.ForwardOptimization = 1; %starts finding the state/parameters backward, -1, (actual time t_k) or Forward, 1 (i.e. at time  t_{k-(w-1)*Nts}

DynOpt.OptXstory = zeros(length(DynOpt.X ),length(time));
DynOpt.OptErrorStory = DynOpt.OptXstory;
DynOpt.OptXstoryTRUE = [state(1:4,:);DynOpt.gamma0*ones(1,length(time));DynOpt.gamma1*ones(1,length(time))];
Jstory = zeros(1,length(time));
J = 1E3;
temp_time = [];

if(DynOpt.TestDynamics ~= 0 )
    DynOpt.Nts = 1;
    DynOpt.U =  zeros(2,2); %  we need to maintain the previous time input value for integration.
    if(DynOpt.TestDynamics == 1) %forward simulation
        DynOpt.X = [z(1);DynOpt.gamma0;DynOpt.gamma1];
        DynOpt.OptXstory(:,1) = DynOpt.X;
        DynOpt.U(:,end) = [int_z(1);elong(1)];%[Dato.ALHcalc;Dato.elong];
        for k=2:length(time),
            %updating the moving window data with the inter-sampled data
            DynOpt.U(:,1:end-1) = DynOpt.U(:,2:end);
            DynOpt.U(:,end) = [int_z(k);elong(k)];%[Dato.ALHcalc;Dato.elong]; 
            DynOpt.X = PlantJumpMap(DynOpt.X,1,1,1);
            DynOpt.OptXstory(:,k) = DynOpt.X;
        end
    else %backward in time
        DynOpt.X = [state(1:4,end);gamma0;gamma1];
        DynOpt.OptXstory(:,length(time)) = DynOpt.X;
        DynOpt.U = [int_z(:)';elong(:)'];%[Dato.ALHcalc;Dato.elong];
        DynOpt.WindowSamples = length(DynOpt.U(1,:));
        for k=1:length(time)-1,
            DynOpt.X = PlantJumpMap(DynOpt.X,DynOpt.WindowSamples-k,-1,1);
            DynOpt.OptXstory(:,length(time)-k) = DynOpt.X;
        end
    
    end
   
else
    
    disp('Processing data with the optimization-based observer...')
    tic
    for k=1:length(time),
        %updating the moving window data with the inter-sampled data
        DynOpt.U(:,1:end-1) = DynOpt.U(:,2:end);
        DynOpt.U(:,end) = [int_z(k);elong(k)];%[Dato.ALHcalc;Dato.elong];
        DynOpt.ActualTimeIndex = k;%it might be useful for non-stationary plants
        DynOpt.Xtrue = [state(:,k);DynOpt.gamma0;DynOpt.gamma1];
        if(k>1)%forward propagation of the previous estimate
            DynOpt.X = PlantJumpMap(DynOpt.X,DynOpt.WindowSamples-1,1,1);
            DynOpt.OptXstory(:,k) = DynOpt.X;
        end
        
        if(mod(k,DynOpt.Nts)-1 == 0)
            % Display iteration step
            disp(['Iteration Number: ', num2str(floor(k/(DynOpt.Nts+1))),'/',num2str(floor(length(time)/(DynOpt.Nts+1)))])
            
            %%%% OUTPUT measurements
            DynOpt.Y(:,1:end-1) = DynOpt.Y(:,2:end);
            DynOpt.Y(:,end) = EvaluateCostFunctionOnWindow_Output_v1(state(:,k),k,1);
            temp_time = [temp_time k];
            %%%%
            
            if( k >= max(1,DynOpt.WindowSamples) )%enough samples have been acquired
                
                if(DynOpt.ForwardOptimization == 1) %forward optimization
                    BackTimeIndex = k-(DynOpt.w-1)*DynOpt.Nts; 
                    [NewXopt, J] = fminunc(@EvaluateCostFunctionOnWindow_v1,DynOpt.OptXstory(:,BackTimeIndex),myoptioptions);
                    %[NewXopt, J] = fmincon(@EvaluateCostFunctionOnWindow_v1,DynOpt.OptXstory(:,BackTimeIndex),DynOpt.Acon,DynOpt.Bcon,DynOpt.Aeq,DynOpt.Beq,DynOpt.lb,DynOpt.ub,@mycon,myoptioptions);
                    DynOpt.X = NewXopt;
                    DynOpt.OptXstory(:,BackTimeIndex) = DynOpt.X;
                    for j =1:DynOpt.WindowSamples-1,%in absolute time as: k-(w-1)*Nts+1:k,
                        DynOpt.X =  PlantJumpMap(DynOpt.X,j,1,1);
                        DynOpt.OptXstory(:,BackTimeIndex + j) = DynOpt.X;
                     end
                    
                else %backward optimization
                    [NewXopt, J] = fminunc(@EvaluateCostFunctionOnWindow_v1,DynOpt.OptXstory(:,k),myoptioptions);
                    DynOpt.X = NewXopt;
                end
            end
            
        end
        Jstory(k) = J;
        DynOpt.OptXstory(:,k) = DynOpt.X;
        clc;
        DynOpt.OptErrorStory(:,k) = DynOpt.Xtrue - DynOpt.X;
    end
    toc
end
%%

figure('Name','J index')
semilogy(time,Jstory)

figure('Name',['Estimates shots '])
for k=1:length(DynOpt.X),
    ax(k)=subplot(length(DynOpt.X),1,k);
    eval(['plot(time,DynOpt.OptXstory(' num2str(k) ',:),''-'',time,DynOpt.OptXstoryTRUE(' num2str(k) ',:),'':'',''LineWidth'' ,2)']);
    eval(['ylabel(''X_' num2str(k) ''')' ]);
    grid on
end
linkaxes(ax,'x');

figure('Name','Estimation Errors')
for k=1:length(DynOpt.X),
    ax(k)=subplot(length(DynOpt.X),1,k);
    eval(['plot(time,DynOpt.OptXstoryTRUE(' num2str(k) ',:)-DynOpt.OptXstory(' num2str(k) ',:),''LineWidth'' ,2)']);
    eval(['ylabel(''E_' num2str(k) ''')' ]);
    grid on
end
linkaxes(ax,'x');

%%
figure('Name','Check windowed data')
clear ax
WindowTime = time(temp_time(end-DynOpt.w+1:1:end))
ax(1) = subplot(2,1,1)
plot(WindowTime,DynOpt.Y(1,:),'s',time,state(2,:),'--','LineWidth',2,'MarkerSize',5)


