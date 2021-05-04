%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STORING TEST DATA

i = endOfTimeSimulation;

%Store true and estimated data for all experiments (test)
ObserverTest.AllSimulation(ObserverTest.Ntest).iner_ECI = ObserverTest.satellites_iner_ECI_alltimes;
ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_iner_ECI = ObserverTest.estimated_satellites_iner_ECI;
ObserverTest.AllSimulation(ObserverTest.Ntest).attitude = satellites_attitude_out;
ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_attitude = ObserverTest.estimated_satellites_attitude;


window_interval = max(1,floor(ObserverTest.StartIntervalWindowPercentage*(i))):1:floor(ObserverTest.EndIntervalWindowPercentage*(i));
ObserverTest.window(Ntest,:) = time([window_interval(1),window_interval(end)]);

for k=1:ObserverTest.Nagents,
    ObserverTest.Mean(Ntest,k) = mean(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.Sigma(Ntest,k) =  std(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.maxError(Ntest,k) =  max(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.TransientTime(Ntest,k) = find(Agent(k).NormEstimationErrorXYZ <= ObserverTest.Mean(Ntest,k),1);
    ObserverTest.MeanDistances(Ntest,k) = mean(Agent(k).ErrorAposterioriDistances(window_interval));
    ObserverTest.SigmaDistances(Ntest,k) =  std(Agent(k).ErrorAposterioriDistances(window_interval));

end

if(ObserverTest.SaveData)
    if(ObserverTest.CentralizedOptimization==1)
        save(['Observer/CentralizedTest/' ObserverTest.FileName],'ObserverTest');
    else
        save(['Observer/DecentralizedTest/' ObserverTest.FileName],'ObserverTest');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


Nfigure = 20;
MyFontSize = 18;
nF = 1.15;
satn  = 3;
ttime = [time(1):time_step:time(i)];


% dx = mypseudo_derivative( ttime, satellites_iner_ECI_alltimes(6*(satn-1)+1,1:i),2,6);
% dy = mypseudo_derivative( ttime, satellites_iner_ECI_alltimes(6*(satn-1)+2,1:i),2,6);
% dz = mypseudo_derivative( ttime, satellites_iner_ECI_alltimes(6*(satn-1)+3,1:i),2,6);
% dxhat = mypseudo_derivative( ttime, xHatUKF(6*(satn-1)+1,1:i),2,6);
% dyhat = mypseudo_derivative( ttime, xHatUKF(6*(satn-1)+2,1:i),2,6);
% dzhat = mypseudo_derivative( ttime, xHatUKF(6*(satn-1)+3,1:i),2,6);
ksatn=satn;
mytime = params.time_step*[1:1:i];

QuaternionError = zeros(ObserverTest.Nagents,i);
OmegaError = zeros(ObserverTest.Nagents,i);
MagBiasError = zeros(ObserverTest.Nagents,i);
GyroBiasError = zeros(ObserverTest.Nagents,i);
for k=1:ObserverTest.Nagents,  
    NumericAgent(k).FilteredErrorXYZ = sqrt( (Agent(k).iner_ECI(1,1:i) - Agent(k).GPSfiltered(1,1:i)).^2 + ...
        (Agent(k).iner_ECI(2,1:i) - Agent(k).GPSfiltered(2,1:i)).^2 + ...
        (Agent(k).iner_ECI(3,1:i) - Agent(k).GPSfiltered(3,1:i)).^2);
    NumericAgent(k).FilteredErrorDXDYDZ = sqrt( (Agent(k).iner_ECI(4,1:i) - Agent(k).GPSderivative(1,1:i)).^2 + ...
        (Agent(k).iner_ECI(5,1:i) - Agent(k).GPSderivative(2,1:i)).^2 + ...
        (Agent(k).iner_ECI(6,1:i) - Agent(k).GPSderivative(3,1:i)).^2);
     NumericAgent(k).GPSpeedError = sqrt( (Agent(k).iner_ECI(4,1:i) - Agent(k).GPSpeed(1,1:i)).^2 + ...
        (Agent(k).iner_ECI(5,1:i) - Agent(k).GPSpeed(2,1:i)).^2 + ...
        (Agent(k).iner_ECI(6,1:i) - Agent(k).GPSpeed(3,1:i)).^2);
    if(ObserverTest.AttitudeObserverOn)
        
        for j=1:i
            temp = quat2angle(Agent(k).attitude_error(1:4,j)');
            QuaternionError(k,j) =  norm(temp);
        end
        OmegaError(k,:) = sqrt(Agent(k).attitude_error(5,1:i).^2 + Agent(k).attitude_error(6,1:i).^2+Agent(k).attitude_error(7,1:i).^2);
        MagBiasError(k,:) = sqrt(Agent(k).magyrobias_error(1,1:i).^2+Agent(k).magyrobias_error(2,1:i).^2+Agent(k).magyrobias_error(3,1:i).^2);
        GyroBiasError(k,:) = sqrt(Agent(k).magyrobias_error(4,1:i).^2+Agent(k).magyrobias_error(5,1:i).^2+Agent(k).magyrobias_error(6,1:i).^2);
    end
end
%%

if(ObserverTest.AttitudeObserverOn)
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    ax(1) = subplot(4,1,1);
    semilogy(ttime,abs(QuaternionError),'LineWidth',2);
    ylabel('q error','FontSize',MyFontSize,'Interpreter','Latex');
    set(gca,'FontSize',MyFontSize);
    ax(2) = subplot(4,1,2);
    semilogy(ttime,abs(OmegaError),'LineWidth',2);
    ylabel('$\omega$ error','FontSize',MyFontSize,'Interpreter','Latex');
    set(gca,'FontSize',MyFontSize);
    ax(3) = subplot(4,1,3);
    semilogy(ttime,abs(MagBiasError),'LineWidth',2);
    ylabel('Mag bias error','FontSize',MyFontSize,'Interpreter','Latex');
    set(gca,'FontSize',MyFontSize);
    ax(4) = subplot(4,1,4);
    semilogy(ttime,abs(GyroBiasError),'LineWidth',2);
    ylabel('$\omega$ bias error','FontSize',MyFontSize,'Interpreter','Latex');
    linkaxes(ax,'x');
    set(gca,'FontSize',MyFontSize);
    
    %%
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    Angles = zeros(3,length(mytime));
    Anglesest = zeros(3,length(mytime));
    AngleErrors = zeros(3,length(mytime));
    for j=1:length(mytime),
        if(j==1)
            [r,p,y] = quat2angle(Agent(satn).attitude(1:4,j)');
        else
            [r,p,y] = mysmooth_quat2angle(Agent(satn).attitude(1:4,j)',Angles(:,j-1));
        end
        Angles(:,j) = [r;p;y];
        if(j==1)
            [re,pe,ye] = quat2angle(Agent(satn).attitude_xHatUKF(1:4,j)');
        else
            [re,pe,ye] = mysmooth_quat2angle(Agent(satn).attitude_xHatUKF(1:4,j)',Anglesest(:,j-1));
        end
        Anglesest(:,j) = [re;pe;ye];
        zeros(3,length(mytime));
        AngleErrors(:,j) = Angles(:,j) - Anglesest(:,j);
    end
    ax(1) = subplot(3,1,1);
    plot(mytime,Angles(1,:)*180/pi,'-',mytime,Anglesest(1,:)*180/pi,'--','LineWidth',2);
    ylabel('Roll','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    title(['Agent '  num2str(satn) ': solid true, dashed estimated'])
    ax(2) = subplot(3,1,2);
    plot(mytime,Angles(2,:)*180/pi,'-',mytime,Anglesest(2,:)*180/pi,'--','LineWidth',2);
    ylabel('Pitch','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    ax(3) = subplot(3,1,3);
    plot(mytime,Angles(3,:)*180/pi,'-',mytime,Anglesest(3,:)*180/pi,'--','LineWidth',2);
    ylabel('Yaw','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    linkaxes(ax,'x');
    set(gca,'FontSize',MyFontSize);
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    ax(1) = subplot(3,1,1);
    plot(mytime,AngleErrors(1,:)*180/pi,'-','LineWidth',2);
    ylabel('Roll Error','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    title(['Agent '  num2str(satn) ': Euler angle estimation error'])
    ax(2) = subplot(3,1,2);
    plot(mytime,AngleErrors(2,:)*180/pi,'-','LineWidth',2);
    ylabel('Pitch Error','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    ax(3) = subplot(3,1,3);
    plot(mytime,AngleErrors(3,:)*180/pi,'-','LineWidth',2);
    ylabel('Yaw Error','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    linkaxes(ax,'x');
    set(gca,'FontSize',MyFontSize);
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    ax(1) = subplot(2,1,1);
    plot(mytime,Agent(satn).magyrobias_error(1:3,1:i),'-','LineWidth',2);
    ylabel('MAG bias error','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    title(['Agent '  num2str(satn) ': Bias errors'])
    ax(2) = subplot(2,1,2);
    plot(mytime,Agent(satn).magyrobias_error(4:end,1:i),'-','LineWidth',2);
    ylabel('Gyro bias rror','FontSize',MyFontSize,'Interpreter','Latex');
    grid on
    linkaxes(ax,'x');
    set(gca,'FontSize',MyFontSize);
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    plot(mytime,Agent(satn).magneto(1:3,1:i),'-','LineWidth',2);
    grid on
    title(['Agent '  num2str(satn) ': Sensed Magnetic Field Vector'])
    set(gca,'FontSize',MyFontSize);
    
    
    
    
    %%
end
%%



%%
figure(Nfigure)
Nfigure = Nfigure + 1;

Angles = zeros(3,length(mytime));
for j=1:length(mytime),
    if(j==1)
            [r,p,y] = quat2angle(Agent(satn).attitude(1:4,j)');
    else
            [r,p,y] = mysmooth_quat2angle(Agent(satn).attitude(1:4,j)',Angles(:,j-1));
    end
    Angles(:,j) = [r;p;y];
end
ax(1) = subplot(3,1,1);
plot(mytime,Angles(1,:)*180/pi,'-','LineWidth',2);
ylabel('Roll','FontSize',MyFontSize,'Interpreter','Latex');
grid on
title(['Agent '  num2str(satn) ': estimated Euler angles'])
ax(2) = subplot(3,1,2);
plot(mytime,Angles(2,:)*180/pi,'LineWidth',2);
ylabel('Pitch','FontSize',MyFontSize,'Interpreter','Latex');
grid on
ax(3) = subplot(3,1,3);
plot(mytime,Angles(3,:)*180/pi,'LineWidth',2);
ylabel('Yaw','FontSize',MyFontSize,'Interpreter','Latex');
grid on
linkaxes(ax,'x');
set(gca,'FontSize',MyFontSize);




figure(Nfigure)
Nfigure = Nfigure + 1;
for j=1:3,
    ax(j) = subplot(3,1,j);
    plot(mytime,Agent(1).iner_ECI(j,1:i),'-k',mytime,Agent(1).GPS(j,1:i),'r-',mytime,Agent(1).GPSopt(j,1:i),'b-','LineWidth',2);
    if(j==1)
        title('GPS on single agent')
    end
    ylim([0.7*min(Agent(1).iner_ECI(j,1:i)), 1.2*max(Agent(1).iner_ECI(j,1:i))]);
    grid on
end
 legend('true','GPS','GPSopt')
set(gca,'FontSize',MyFontSize);



figure(Nfigure)
Nfigure = Nfigure + 1;
for j=1:3,
    ax(j) = subplot(3,1,j);
    plot(mytime,Agent(1).iner_ECI(j,1:i),'-k',mytime,Agent(1).GPSfiltered(j,1:i),'-g',mytime,Agent(1).GPS(j,1:i),'r-','LineWidth',2);
    if(j==1)
        title('Filtered position')
    end
    grid on
end
legend('True','Filtered','GPS')
set(gca,'FontSize',MyFontSize);

figure(Nfigure)
Nfigure = Nfigure + 1;
for j=1:3,
    ax(j) = subplot(3,1,j);
    plot(mytime,Agent(k).GPSderivative(j,1:i),'--',mytime,Agent(k).GPSpeed(j,1:i),':',mytime,Agent(k).iner_ECI(j+3,1:i),'-','LineWidth',2);
    if(j==1)
        title('Speed estimation ')
    end
    grid on
end
legend('Used','GPSpeed','True')
set(gca,'FontSize',MyFontSize);


figure(Nfigure)
Nfigure = Nfigure + 1;

clear ax
temp = zeros(ObserverTest.Nagents,i);
for k=1:ObserverTest.Nagents,
    temp(k,:) = NumericAgent(k).FilteredErrorXYZ;
end
ax(1) = subplot(3,1,1);
semilogy(ttime,temp,'LineWidth',2);
grid on
ylabel('$||E_i(x,y,z)|| Filter $','FontSize',MyFontSize,'Interpreter','Latex')
set(gca,'FontSize',MyFontSize);

ax(2) = subplot(3,1,2);
temp = zeros(ObserverTest.Nagents,i);
for k=1:ObserverTest.Nagents,
    temp(k,:) = NumericAgent(k).FilteredErrorDXDYDZ;
end
semilogy(ttime,temp,'LineWidth',2);
grid on
ylabel('$||dE_i(x,y,z)||, \hat{v} $','FontSize',MyFontSize,'Interpreter','Latex')
set(gca,'FontSize',MyFontSize);

ax(3) = subplot(3,1,3);
temp = zeros(ObserverTest.Nagents,i);
for k=1:ObserverTest.Nagents,
    temp(k,:) = NumericAgent(k).GPSpeedError;
end
semilogy(ttime,temp,'LineWidth',2);
grid on
ylabel('$||dE_i(x,y,z)||, v_{GPS} $','FontSize',MyFontSize,'Interpreter','Latex')
xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
set(gca,'FontSize',MyFontSize);


% decrease distance between subplots
for h = 1:length(ax)
    vCurrPos = get(ax(h), 'position'); % current position
    %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
    set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
    if h < length(ax)
        set(ax(h), 'XTickLabel', ' ')
    end
end
linkaxes(ax,'x');


if(ObserverTest.ObserverON)
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    ax(1) = subplot(3,1,1)
    plot(ttime,Agent(ksatn).xHatUKF(1,1:i)-Agent(1).iner_ECI(1,1)*ones(1,length(ttime)),'r-',...
        ttime,Agent(ksatn).xHatUKF(2,1:i),'b-',...
        ttime, Agent(ksatn).xHatUKF(3,1:i),'c-','LineWidth',2);
    legend('x','y','z','FontSize',MyFontSize,'Interpreter','Latex')
    hold on
    %plot(ttime, myEulerCopy(6*(satn-1)+1:6*(satn-1)+3,1:i)','--','LineWidth',2);
    plot(ttime, Agent(ksatn).iner_ECI(1,1:i)-Agent(1).iner_ECI(1,1)*ones(1,length(ttime)),'r--',...
        ttime,Agent(ksatn).iner_ECI(2,1:i)','b--',...
        ttime,Agent(ksatn).iner_ECI(3,1:i)','c--','LineWidth',2);
    hold off
    grid on
    ylabel('[km]')
    
    ax(2) = subplot(3,1,2)
    plot(ttime,Agent(ksatn).xHatUKF(4,1:i),'r-',...
        ttime,Agent(ksatn).xHatUKF(5,1:i),'b-',...
        ttime,Agent(ksatn).xHatUKF(6,1:i),'c-','LineWidth',2);
    legend('$\dot x$','$\dot y$','$\dot z$','FontSize',MyFontSize,'Interpreter','Latex')
    hold on
    set(gca,'FontSize',MyFontSize);
    %plot(ttime, myEulerCopy(6*(satn-1)+4:6*(satn-1)+6,1:i)','--','LineWidth',2);
    plot(ttime, Agent(ksatn).iner_ECI(4,1:i)','r--',...
        ttime,  Agent(ksatn).iner_ECI(5,1:i)','b--',...
        ttime,  Agent(ksatn).iner_ECI(6,1:i)','c--','LineWidth',2);
    hold off
    grid on
    xlim([time(1) time(i)])
    xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
    ylabel('[km/s]')
    set(gca,'FontSize',MyFontSize);
    
    ax(3) = subplot(3,1,3);
    plot(ttime,Agent(ksatn).iner_ECI(1,1:i),'r-',...
        ttime,Agent(ksatn).iner_ECI(2,1:i),'b-',...
        ttime, Agent(ksatn).iner_ECI(3,1:i),'c-','LineWidth',2);
    legend('dx','dy','dz','FontSize',MyFontSize,'Interpreter','Latex')
    grid on
    ylabel('[m]')
    
   
    % decrease distance between subplots
    for h = 1:length(ax)
        vCurrPos = get(ax(h), 'position'); % current position
        %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
        set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
        if h < length(ax)
            set(ax(h), 'XTickLabel', ' ')
        end
    end
    linkaxes(ax,'x');
    
    
    %%
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    ax(1) = subplot(2,1,1)
    plot(ttime,Agent(ksatn).iner_ECI_EstimationError(1,1:i),'r-',...
        ttime,Agent(ksatn).iner_ECI_EstimationError(2,1:i),'b-',...
        ttime,Agent(ksatn).iner_ECI_EstimationError(3,1:i),'c-','LineWidth',2);
    grid on
    ylabel('[m]')
    legend('$E_x$','$E_y$','$E_z$','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    ax(2) = subplot(2,1,2)
    plot(ttime,Agent(ksatn).iner_ECI_EstimationError(4,1:i),'r-',...
        ttime,Agent(ksatn).iner_ECI_EstimationError(5,1:i),'b-',...
        ttime,Agent(ksatn).iner_ECI_EstimationError(6,1:i),'c-','LineWidth',2);
    grid on
    xlim([time(1) time(i)])
    legend('$E_{\dot x}$','$E_{\dot y}$','$E_{\dot z}$','FontSize',MyFontSize,'Interpreter','Latex')
    ylabel('[m/s]')
    set(gca,'FontSize',MyFontSize);
    % decrease distance between subplots
    for h = 1:length(ax)
        vCurrPos = get(ax(h), 'position'); % current position
        %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
        set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
        if h < length(ax)
            set(ax(h), 'XTickLabel', ' ')
        end
    end
    linkaxes(ax,'x');
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    temp = zeros(ObserverTest.Nagents,i);
    for ksatn=1:ObserverTest.Nagents,
        temp(ksatn,:) = Agent(ksatn).NormEstimationErrorXYZ(1:i);
    end
    
    clear ax
    ax(1) = subplot(3,1,1);
    semilogy(ttime,temp,'LineWidth',2);
    hold on
    semilogy(ttime, norm(ObserverTest.GPSGaussianCovariance(1:3))*ones(1,i),'k--','LineWidth',2);
    hold off
    grid on
    ylabel('$||E_i(x,y,z)||$','FontSize',MyFontSize,'Interpreter','Latex')
    
    ax(2) = subplot(3,1,2);
    temp = zeros(ObserverTest.Nagents,i);
    for ksatn=1:ObserverTest.Nagents,
        temp(ksatn,:) = sqrt(Agent(ksatn).iner_ECI_EstimationError(4,1:i).^2+Agent(ksatn).iner_ECI_EstimationError(5,1:i).^2+...
        Agent(ksatn).iner_ECI_EstimationError(6,1:i).^2);
        %temp(ksatn,1:i) = Agent(ksatn).UWBcorrectionActive(1:i)';
    end
    semilogy(ttime,temp,'LineWidth',2)
    grid on
    ylabel('$||dE_i(dx,dy,dz)||$','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    
    
    ax(3) = subplot(3,1,3);
    temp = zeros(ObserverTest.Nagents,i);
    for ksatn=1:ObserverTest.Nagents,
        temp(ksatn,1:i) = Agent(ksatn).MismatchAprioriDistances(1:i)';
    end
    plot(ttime,mean(temp),'-b','LineWidth',2)
    hold on
    temp = zeros(ObserverTest.Nagents,i);
    for ksatn=1:ObserverTest.Nagents,
        temp(ksatn,1:i) = Agent(ksatn).ErrorAprioriDistances(1:i)';
    end
    plot(ttime,mean(temp),':b','LineWidth',2)
    temp = zeros(ObserverTest.Nagents,i);
    for ksatn=1:ObserverTest.Nagents,
        temp(ksatn,1:i) = Agent(ksatn).MismatchAposterioriDistances(1:i)';
    end
    plot(ttime,mean(temp),'-r','LineWidth',3)
    
    temp = zeros(ObserverTest.Nagents,i);
    for ksatn=1:ObserverTest.Nagents,
        temp(ksatn,1:i) = Agent(ksatn).ErrorAposterioriDistances(1:i)';
    end
    plot(ttime,mean(temp),':r','LineWidth',3)
    hold off
    grid on
    %ylim([0 10])
    ylabel('SumDistances','FontSize',MyFontSize,'Interpreter','Latex')
    legend('Pre mismatch','Pre error','Post mismatch','Post error','FontSize',MyFontSize);
    xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    
    % decrease distance between subplots
    for h = 1:length(ax)
        vCurrPos = get(ax(h), 'position'); % current position
        %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
        set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
        if h < length(ax)
            set(ax(h), 'XTickLabel', ' ')
        end
    end
    linkaxes(ax,'x');
    set(gca,'FontSize',MyFontSize);
    
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.Mean(Ntest,:)*1000,ObserverTest.Sigma(Ntest,:)*1000,'LineWidth',2);
    ylabel('Position: mean+-sigma [m]','FontSize',MyFontSize,'Interpreter','Latex')
    xlabel('Agents','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    grid on
    
    figure(Nfigure)
    Nfigure = Nfigure + 1;
    
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.MeanDistances(Ntest,:)*1000,ObserverTest.SigmaDistances(Ntest,:)*1000,'b--','LineWidth',2);
    ylabel('Relative Distances: mean+-sigma [m]','FontSize',MyFontSize,'Interpreter','Latex')
    xlabel('Agents','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    grid on
    
    if(ObserverTest.UWBoptimizationOn == 1)
      %%  
        figure(Nfigure)
        Nfigure = Nfigure + 1;
        
        temp1 = zeros(ObserverTest.Nagents,i);
        temp2 = zeros(ObserverTest.Nagents,i);
        temp3 = zeros(ObserverTest.Nagents,i);
        temp4 = zeros(ObserverTest.Nagents,i);
        temp5 = zeros(ObserverTest.Nagents,i);
        for ksatn=1:ObserverTest.Nagents,
            temp = Agent(ksatn).iner_ECI([1 2 3],1:i) - Agent(ksatn).GPS(:,1:i);
            temp1(ksatn,:) = sqrt(sum(temp.*temp));
            temp = Agent(ksatn).iner_ECI([1 2 3],1:i) - Agent(ksatn).GPSopt(:,1:i);
            temp2(ksatn,:) = sqrt(sum(temp.*temp));   
            temp = Agent(ksatn).iner_ECI([4 5 6],1:i) - Agent(ksatn).GPSpeed(:,1:i);
            temp3(ksatn,:) = sqrt(sum(temp.*temp));            
            temp = Agent(ksatn).iner_ECI([4 5 6],1:i) - Agent(ksatn).GPSpeedOpt(:,1:i);
            temp4(ksatn,:) = sqrt(sum(temp.*temp));      
            temp = Agent(ksatn).GPS(:,1:i) - Agent(ksatn).GPSopt(:,1:i);      
            temp5(ksatn,:) = sqrt(sum(temp.*temp));
        end
        clear ax
        ax(1) = subplot(5,1,1);
        plot(ttime,temp1,'LineWidth',2);
        hold off
        grid on
        ylabel('E_{GPS}','FontSize',MyFontSize)
        
        ax(2) = subplot(5,1,2);
        plot(ttime,temp2,'LineWidth',2);
        hold off
        grid on
        ylabel('E_{GPS,opt} ','FontSize',MyFontSize)
        
        ax(3) = subplot(5,1,3);
        semilogy(ttime, sum(temp1),ttime, sum(temp2),'LineWidth',2);
        hold off
        grid on
        legend('$\Sigma ||E_{i,GPS}||$','$\Sigma|| E_{i,GPSopt}||$','FontSize',MyFontSize)
        set(gca,'FontSize',MyFontSize);
        
        ax(4) = subplot(5,1,4);
        plot(ttime, temp5,'LineWidth',2);
        hold off
        grid on
        legend('$GPS-GPS_{opt}$','FontSize',MyFontSize)
        set(gca,'FontSize',MyFontSize);
        
        
        ax(5) = subplot(5,1,5);
        semilogy(ttime, sum(temp3),ttime, sum(temp4),'LineWidth',2);
        hold off
        grid on
        legend('$\Sigma ||E_{i,GPSpeed}||$','$\Sigma||E_{i,GPSpeedopt}||$','FontSize',MyFontSize)
        
        xlabel('Time (s)','FontSize',MyFontSize)
        set(gca,'FontSize',MyFontSize);
        
        % decrease distance between subplots
        for h = 1:length(ax)
            vCurrPos = get(ax(h), 'position'); % current position
            %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
            set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
            if h < length(ax)
                set(ax(h), 'XTickLabel', ' ')
            end
        end
        linkaxes(ax,'x');
        set(gca,'FontSize',MyFontSize);
        %%
        
        
%         figure(Nfigure)
%         Nfigure = Nfigure + 1;
%         set(gca,'FontSize',MyFontSize);
%         %ax(1) = subplot(2,1,1);
%         errorbar(1:length(Test.Mean(1,:)), Test.Mean(Ntest,:),Test.Sigma(Ntest,:),'LineWidth',2);
%         grid on
%         xlabel('Agents','FontSize',MyFontSize,'Interpreter','Latex')
%         ylabel('Mean+/-Sigma')


        
        
    
    
    % % temp = zeros(ObserverTest.Nagents,i);
    % % for k = 1:i,
    % %     for j = 1:ObserverTest.Nagents,
    % %        temp(j,k) = norm(XhatEuler(ExtractGPS+(j-1)*6,k)-satellites_iner_ECI_alltimes(ExtractGPS+(j-1)*6,k));
    % %     end
    % % end
    % % figure(4)
    % % semilogy(ttime,temp,'LineWidth',2);
    % % grid on
    % % ylabel('Euler approximation ||E_i(x,y,z)||','FontSize',MyFontSize)
    % % xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
    % % set(gca,'FontSize',MyFontSize);
    % %%
    % figure(5)
    % plot(ttime,XhatEuler(ExtractGPS+(satn-1)*6,1:i),'-',ttime,satellites_iner_ECI_alltimes(ExtractGPS+(satn-1)*6,1:i),'LineWidth',2);
    % grid on
    % ylabel('E_?(x,y,z)','FontSize',MyFontSize)
    % xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
    % set(gca,'FontSize',MyFontSize);
    %%
    
%     if(Ntest == Test.TotalRuns)
%         figure(Nfigure)
%         Nfigure = Nfigure + 1;
%         colore ={'r','b','k','y','g','m','c','b','y'}
%         
%         for kt = Test.UWBOptimizationNoBeforeThan-3:i ,
%             %         for ksatn=1:ObserverTest.Nagents,
%             %             coldraw = [colore{ksatn} 'x'];
%             %             temp = Agent(ksatn).iner_ECI([1 2 3],kt);
%             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
%             %             hold on
%             %             coldraw = [colore{ksatn} 'o'];
%             %             temp = Agent(ksatn).xHatUKF([1 2 3],kt);
%             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
%             %             coldraw = [colore{ksatn} 's'];
%             %             temp = Agent(ksatn).GPS(:,kt);
%             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
%             %             coldraw = [colore{ksatn} 'd'];
%             %             temp = Agent(ksatn).GPSopt(:,kt);
%             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
%             %         end
%             %         hold off
%             %         grid on
%             %         pause
%             for ksatn=1:ObserverTest.Nagents,
%                 ax(1) = subplot(2,1,1);
%                 coldraw = [colore{ksatn} 'x'];
%                 temp = Agent(ksatn).iner_ECI([1 2 3],kt);
%                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
%                 hold on
%                 temp = Agent(ksatn).iner_ECI([1 2 3],kt-1:kt);
%                 plot(temp(1,:),temp(2,:),[colore{ksatn} '-'],'LineWidth',2);
%                 
%                 coldraw = [colore{ksatn} 'o'];
%                 temp = Agent(ksatn).xHatUKF([1 2 3],kt);
%                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
%                 temp = Agent(ksatn).xHatUKF([1 2 3],kt-1:kt);
%                 plot(temp(1,:),temp(2,:),[colore{ksatn} '--'],'LineWidth',2);
%                 
%                 coldraw = [colore{ksatn} 's'];
%                 temp = Agent(ksatn).GPS(:,kt);
%                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
%                 
%                 coldraw = [colore{ksatn} 'd'];
%                 temp = Agent(ksatn).GPSopt(:,kt);
%                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
%                 
%                 
%                 ax(2) = subplot(2,1,2);
%                 coldraw = [colore{ksatn} 'x'];
%                 temp = Agent(ksatn).iner_ECI([1 2 3],kt);
%                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
%                 hold on
%                 temp = Agent(ksatn).iner_ECI([1 2 3],kt-1:kt);
%                 plot(temp(1,:),temp(3,:),[colore{ksatn} '-'],'LineWidth',2);
%                 
%                 coldraw = [colore{ksatn} 'o'];
%                 temp = Agent(ksatn).xHatUKF([1 2 3],kt);
%                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
%                 temp = Agent(ksatn).xHatUKF([1 2 3],kt-1:kt);
%                 plot(temp(1,:),temp(3,:),[colore{ksatn} '--'],'LineWidth',2);
%                 
%                 
%                 coldraw = [colore{ksatn} 's'];
%                 temp = Agent(ksatn).GPS(:,kt);
%                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
%                 coldraw = [colore{ksatn} 'd'];
%                 temp = Agent(ksatn).GPSopt(:,kt);
%                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
%             end
%             ax(1) = subplot(2,1,1);
%             hold off
%             ax(2) = subplot(2,1,2);
%             hold off
%             if(kt == Test.UWBOptimizationNoBeforeThan-3)
%                 ax(1) = subplot(2,1,1);
%                 xlabel('X','FontSize',MyFontSize);
%                 ylabel('Y','FontSize',MyFontSize);
%                 grid on
%                 title('x true, o estimate, square GPS, diamond GPS opt')
%                 ax(2) = subplot(2,1,2);
%                 xlabel('X','FontSize',MyFontSize);
%                 ylabel('Z','FontSize',MyFontSize);
%                 grid on
%             end
%             
%             if(kt== Test.UWBOptimizationNoBeforeThan-3)
%                 %pause;
%             else
%                 %pause(0.3);
%             end
%         end
%         grid on
%         hold off
%         
%         linkaxes(ax,'x');
%         %%
%     end
    end
    
    
else
    
    figure(Nfigure)
        Nfigure = Nfigure + 1;
        
    ax(1) = subplot(2,1,1)
    plot(ttime,satellites_iner_ECI_alltimes((ksatn-1)*6+1,1:i),'r-',...
        ttime,satellites_iner_ECI_alltimes((ksatn-1)*6+2,1:i),'b-',...
        ttime, satellites_iner_ECI_alltimes((ksatn-1)*6+3,1:i),'c-','LineWidth',2);
    legend('x','y','z','FontSize',MyFontSize,'Interpreter','Latex')
    grid on
    ylabel('[m] (x-x_0)')
    ax(2) = subplot(2,1,2)
    plot(ttime, satellites_iner_ECI_alltimes((ksatn-1)*6+4,1:i),'r--',...
        ttime,satellites_iner_ECI_alltimes((ksatn-1)*6+5,1:i),'b--',...
        ttime,satellites_iner_ECI_alltimes((ksatn-1)*6+6,1:i),'c--','LineWidth',2);
    legend('$\dot x$','$\dot y$','$\dot z$','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    grid on
    xlim([time(1) time(i)])
    xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
    ylabel('[m/s]')
    
    
end
% 
% 
% %% Plots
% % kk = 0;
% % Mean_surf = zeros(Test.N_P,Test.N_R);
% % for np=1:Test.N_P,
% %     for nr=1:Test.N_R,
% %         kk = kk+1;
% %         if(min(Test.UnfeasableR(kk,:))<0.1)
% %             Mean_surf(np,nr) = mean(Test.Mean(kk,:));
% %         else
% %             Mean_surf(np,nr) = 0;
% %         end
% %     end
% % end
% % figure(Nfigure)
% %         Nfigure = Nfigure + 1;
% %         
% % surf(Pindexes,Rindexes,Mean_surf);
% % xlabel('\rho*P')
% % ylabel('\rho*R')
% % grid on

%% 
% figure(Nfigure)
% Nfigure = Nfigure + 1;
% kk = 0;
% minmin = 1E3;
% minmin_kk = 0;
% for j=1:length(Test.N_Wuwb_range),
%     Mean_surf = zeros(length(Test.N_Wgps_range),length(Test.N_Wsigma_range));
%     for ngps=1:length(Test.N_Wgps_range),
%         for nsigma=1:length(Test.N_Wsigma_range),
%             kk = kk+1;
%             if(min(Test.UnfeasableR(kk,:))<0.1)
%                 Mean_surf(ngps,nsigma) = mean(Test.Mean(kk,:));
%                 if(minmin > Mean_surf(ngps,nsigma)+mean(Test.Sigma(kk,:)) )
%                     minmin = Mean_surf(ngps,nsigma) + mean(Test.Sigma(kk,:));
%                     minmin_kk = kk;
%                     Test.Gains(:,minmin_kk)
%                 end
%             else
%                 Mean_surf(ngps,nsigma) = 0;
%             end
%         end
%     end
%     grid on
%     surf(Test.N_Wsigma_range,Test.N_Wgps_range,Mean_surf);
%     zlabel(['$W_{uwb}$ =  ' num2str(Test.N_Wuwb_range(j))],'FontSize',MyFontSize,'Interpreter','Latex')
%     xlabel(['$W_{sigma}$'],'FontSize',MyFontSize,'Interpreter','Latex');
%     ylabel(['$W_{gps}$'],'FontSize',MyFontSize,'Interpreter','Latex');
%     pause
% end
% %%
% figure(Nfigure)
% Nfigure = Nfigure + 1;
% set(gca,'FontSize',MyFontSize);
% %ax(1) = subplot(2,1,1);
% errorbar(1:length(Test.Mean(1,:)), Test.Mean(minmin_kk,:),Test.Sigma(minmin_kk,:),'LineWidth',2);
% grid on
% xlabel('Agents','FontSize',MyFontSize,'Interpreter','Latex')
% ylabel('Mean+/-Sigma')

% % decrease distance between subplots
% for h = 1:length(ax)
%     vCurrPos = get(ax(h), 'position'); % current position
%     %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
%     set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
%     if h < length(ax)
%         set(ax(h), 'XTickLabel', ' ')
%     end
% end
% linkaxes(ax,'x');
% set(gca,'FontSize',MyFontSize);




