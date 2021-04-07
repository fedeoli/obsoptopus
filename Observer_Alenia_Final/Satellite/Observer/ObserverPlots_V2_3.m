%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STORING AND SHOWING TEST DATA

i = length(time);

%Store true and estimated data for all experiments (test)
ObserverTest.AllSimulation(ObserverTest.Ntest).iner_ECI = ObserverTest.satellites_iner_ECI_alltimes;
ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_iner_ECI = ObserverTest.estimated_satellites_iner_ECI;
ObserverTest.AllSimulation(ObserverTest.Ntest).attitude = satellites_attitude_out;
ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_attitude = ObserverTest.estimated_satellites_attitude;


window_interval = max(1,floor(ObserverTest.StartIntervalWindowPercentage*(i))):1:floor(ObserverTest.EndIntervalWindowPercentage*(i));
ObserverTest.window(ObserverTest.Ntest,:) = time([window_interval(1),window_interval(end)]);

for k=1:ObserverTest.Nagents,
    ObserverTest.Mean(ObserverTest.Ntest,k) = mean(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.Sigma(ObserverTest.Ntest,k) =  std(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.maxError(ObserverTest.Ntest,k) =  max(Agent(k).NormEstimationErrorXYZ(window_interval));
    ObserverTest.TransientTime(ObserverTest.Ntest,k) = find(Agent(k).NormEstimationErrorXYZ <= ObserverTest.Mean(ObserverTest.Ntest,k),1);
    ObserverTest.MeanDistances(ObserverTest.Ntest,k) = mean(Agent(k).ErrorAposterioriDistances(window_interval));
    ObserverTest.SigmaDistances(ObserverTest.Ntest,k) =  std(Agent(k).ErrorAposterioriDistances(window_interval));

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
satn  = 3; %for which satellite the estimation has to be shown
ttime = [time(1):time_step:time(i)];


% dx = mypseudo_derivative( ttime, satellites_iner_ECI_alltimes(6*(satn-1)+1,1:i),2,6);
% dy = mypseudo_derivative( ttime, satellites_iner_ECI_alltimes(6*(satn-1)+2,1:i),2,6);
% dz = mypseudo_derivative( ttime, satellites_iner_ECI_alltimes(6*(satn-1)+3,1:i),2,6);
% dxhat = mypseudo_derivative( ttime, xHatUKF(6*(satn-1)+1,1:i),2,6);
% dyhat = mypseudo_derivative( ttime, xHatUKF(6*(satn-1)+2,1:i),2,6);
% dzhat = mypseudo_derivative( ttime, xHatUKF(6*(satn-1)+3,1:i),2,6);
mytime = params.time_step*[1:1:i];

QuaternionError = zeros(ObserverTest.Nagents,i);
OmegaError = zeros(ObserverTest.Nagents,i);
MagBiasError = zeros(ObserverTest.Nagents,i);
GyroBiasError = zeros(ObserverTest.Nagents,i);
Angle_true = zeros(3*ObserverTest.Nagents,i);
Angle_estimated = zeros(3*ObserverTest.Nagents,i);
Angle_errors = zeros(3*ObserverTest.Nagents,i);
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
end
%%

if(ObserverTest.AttitudeObserverOn)
    
    Angle_true = zeros(3*ObserverTest.Nagents,i);
    Angle_estimated = zeros(3*ObserverTest.Nagents,i);
    for k=1:ObserverTest.Nagents,
        Angle_true(3*(k-1)+1:3*(k-1)+3,1) = quat2angle(Agent(k).attitude(1:4,1)');
        Angle_estimated(3*(k-1)+1:3*(k-1)+3,1) = quat2angle(Agent(k).attitude_xHatUKF(1:4,1)');
        for j =2:i,
            [e1, e2, e3] = quat2angle(Agent(k).attitude(1:4,j)');
            e = [e1;e2;e3];
            for r=1:3,
                if( abs(e(r) - Angle_true(3*(k-1)+r,j-1) ) > pi )
                    Angle_true(3*(k-1)+r,j) = e(r) - sign(e(r) - Angle_true(3*(k-1)+r,j-1))*2*pi;
                else
                    Angle_true(3*(k-1)+r,j) = e(r);
                end
            end
            [e1, e2, e3] = quat2angle(Agent(k).attitude_xHatUKF(1:4,j)');
            e = [e1;e2;e3];
            for r=1:3,
                if( abs(e(r) - Angle_estimated(3*(k-1)+r,j-1))  > pi )
                    Angle_estimated(3*(k-1)+r,j) = e(r) - sign(e(r) - Angle_estimated(3*(k-1)+r,j-1))*2*pi;
                else
                    Angle_estimated(3*(k-1)+r,j) = e(r);
                end
            end
        end
        temp = Agent(k).attitude(5:7,1:i) - Agent(k).attitude_xHatUKF(5:7,1:i) ;
        OmegaError(k,:) = sqrt(temp(1,1:i).^2 + temp(2,1:i).^2+temp(3,1:i).^2);
    end
    Angle_errors = Angle_true - Angle_estimated; 
    
    figure('Name', 'Attitude estimation errors')
 
    ax(1) = subplot(4,1,1);
    semilogy(ttime,abs(QuaternionError),'LineWidth',2);
    ylabel('q error')
    
    ax(2) = subplot(4,1,2);
    semilogy(ttime,abs(OmegaError),'LineWidth',2);
    ylabel('$\omega$ error');
    
    ax(3) = subplot(4,1,3);
    semilogy(ttime,abs(MagBiasError),'LineWidth',2);
    ylabel('Mag bias error')%,'FontSize',MyFontSize,'Interpreter','Latex');
    
    ax(4) = subplot(4,1,4);
    semilogy(ttime,abs(GyroBiasError),'LineWidth',2);
    ylabel('$\omega$ bias error')%,'FontSize',MyFontSize,'Interpreter','Latex');
    linkaxes(ax,'x');
    
    %%
    figure('Name', ['Euler angle estimation errors (Sum of abs Errors among agents)'])
    
    temp = [1:3:3*ObserverTest.Nagents];
    ax(1) = subplot(3,1,1);
    plot(mytime,sum(abs(Angle_errors(temp,:)))*180/pi,'-','LineWidth',2);
    ylabel('Roll ');
    grid on
    
    temp = [2:3:3*ObserverTest.Nagents];
    ax(2) = subplot(3,1,2);
    plot(mytime,sum(abs(Angle_errors(temp,:)))*180/pi,'-','LineWidth',2);
    ylabel('Pitch');
    grid on
    
    temp = [3:3:3*ObserverTest.Nagents];
    ax(3) = subplot(3,1,3);
    plot(mytime,sum(abs(Angle_errors(temp,:)))*180/pi,'-','LineWidth',2);
    ylabel('Yaw');
    grid on
    linkaxes(ax,'x');
    
    %%
    figure('Name', ['Agent '  num2str(satn) ' attitude: solid true, dashed estimated'])
    
    ax(1) = subplot(3,1,1);
    plot(mytime,Angle_true(1+3*(satn-1),:)*180/pi,'-',mytime,Angle_estimated(1+3*(satn-1),:)*180/pi,'--','LineWidth',2);
    ylabel('Roll')
    grid on
    ax(2) = subplot(3,1,2);
    plot(mytime,Angle_true(2+3*(satn-1),:)*180/pi,'-',mytime,Angle_estimated(2+3*(satn-1),:)*180/pi,'--','LineWidth',2);
    ylabel('Pitch');
    grid on
    ax(3) = subplot(3,1,3);
    plot(mytime,Angle_true(3+3*(satn-1),:)*180/pi,'-',mytime,Angle_estimated(3+3*(satn-1),:)*180/pi,'--','LineWidth',2);
    ylabel('Yaw');
    grid on
    linkaxes(ax,'x');

    %%
    
    
    figure('Name', ['Agent '  num2str(satn) ' attitude: Euler angle estimation error'])
    
    ax(1) = subplot(3,1,1);
    plot(mytime,Angle_errors(1+3*(satn-1),:)*180/pi,'-','LineWidth',2);
    ylabel('Roll Error');
    grid on
    
    ax(2) = subplot(3,1,2);
    plot(mytime,Angle_errors(2+3*(satn-1),:)*180/pi,'-','LineWidth',2);
    ylabel('Pitch Error');
    grid on
    ax(3) = subplot(3,1,3);
    plot(mytime,Angle_errors(3+3*(satn-1),:)*180/pi,'-','LineWidth',2);
    ylabel('Yaw Error');
    grid on
    linkaxes(ax,'x');
    
    
    figure('Name', ['Agent '  num2str(satn) ' attitude: Bias errors'])
    
    ax(1) = subplot(2,1,1);
    plot(mytime,Agent(satn).magyrobias_error(1:3,1:i),'-','LineWidth',2);
    ylabel('MAG bias error');
    grid on
    
    ax(2) = subplot(2,1,2);
    plot(mytime,Agent(satn).magyrobias_error(4:end,1:i),'-','LineWidth',2);
    ylabel('Gyro bias rror');
    grid on
    linkaxes(ax,'x');
    
    
    figure('Name', ['Agent '  num2str(satn) ' attitude: Sensed Magnetic Vector Fields (-- 2nd) '])
    
    plot(mytime,Agent(satn).magneto(1:3,1:i),'-',mytime,Agent(satn).magneto2(1:3,1:i),'--','LineWidth',2);
    grid on
    
    %%
    figure('Name', ['Agent '  num2str(satn) ' attitude: Sensed Magnetic Field Vector'])
    
    ax(1) = subplot(3,1,1);
    plot(mytime,Agent(satn).gyros(1,1:i),':',mytime,Agent(satn).attitude(5,1:i),'-',mytime,Agent(satn).attitude_xHatUKF(5,1:i),'--','LineWidth',2);
    ylabel('$\Omega_x$');
    legend('measure','true','estimated')
    grid on
    
    ax(2) = subplot(3,1,2);
    plot(mytime,Agent(satn).gyros(2,1:i),':',mytime,Agent(satn).attitude(6,1:i),'-',mytime,Agent(satn).attitude_xHatUKF(6,1:i),'--','LineWidth',2);
    ylabel('$\Omega_y$');
    legend('measure','true','estimated')
    grid on
    
    ax(2) = subplot(3,1,3);
    plot(mytime,Agent(satn).gyros(3,1:i),':',mytime,Agent(satn).attitude(7,1:i),'-',mytime,Agent(satn).attitude_xHatUKF(7,1:i),'--','LineWidth',2);
    ylabel('$\Omega_z$');
    legend('measure','true','estimated')
    grid on
    linkaxes(ax,'x');
    
    
    
    
      %%    
end
% 
% figure('Name', ['Agent '  num2str(satn) ' GPS (opt)'])
% 
% for j=1:3,
%     ax(j) = subplot(3,1,j);
%     plot(mytime,Agent(satn).iner_ECI(j,1:i),'-k',mytime,Agent(satn).GPS(j,1:i),'r-',mytime,Agent(satn).GPSopt(j,1:i),'b-','LineWidth',2);
%     ylim([0.7*min(Agent(satn).iner_ECI(j,1:i)), 1.2*max(Agent(satn).iner_ECI(j,1:i))]);
%     grid on
% end
% legend('true','GPS','GPSopt')
% 
% 
% figure('Name', ['Agent '  num2str(satn) ' velocities '])
% 
% for j=1:3,
%     ax(j) = subplot(3,1,j);
%     plot(mytime,Agent(satn).GPSderivative(j,1:i),'--',mytime,Agent(satn).GPSpeed(j,1:i),':',mytime,Agent(satn).iner_ECI(j+3,1:i),'-','LineWidth',2);
%     grid on
% end
% legend('Used','GPSpeed','True')
% 
% 
% figure('Name', ['Agent '  num2str(satn) ' velocities '])
% 
% clear ax
% temp = zeros(ObserverTest.Nagents,i);
% for k=1:ObserverTest.Nagents,
%     temp(k,:) = NumericAgent(k).FilteredErrorXYZ;
% end
% ax(1) = subplot(3,1,1);
% semilogy(ttime,temp,'LineWidth',2);
% grid on
% ylabel('$||E_i(x,y,z)|| Filter $')
% 
% ax(2) = subplot(3,1,2);
% temp = zeros(ObserverTest.Nagents,i);
% for k=1:ObserverTest.Nagents,
%     temp(k,:) = NumericAgent(k).FilteredErrorDXDYDZ;
% end
% semilogy(ttime,temp,'LineWidth',2);
% grid on
% ylabel('$||dE_i(x,y,z)||, \hat{v} $','FontSize',MyFontSize,'Interpreter','Latex')
% set(gca,'FontSize',MyFontSize);
% 
% ax(3) = subplot(3,1,3);
% temp = zeros(ObserverTest.Nagents,i);
% for k=1:ObserverTest.Nagents,
%     temp(k,:) = NumericAgent(k).GPSpeedError;
% end
% semilogy(ttime,temp,'LineWidth',2);
% grid on
% ylabel('$||dE_i(x,y,z)||, v_{GPS} $','FontSize',MyFontSize,'Interpreter','Latex')
% xlabel('Time (s)','FontSize',MyFontSize,'Interpreter','Latex')
% set(gca,'FontSize',MyFontSize);
% 
% 
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
% 

if(ObserverTest.ObserverON)
    
   figure('Name', ['Agent '  num2str(satn) ': state (true--, estimated -) '])
   
    ax(1) = subplot(2,1,1);
    plot(ttime,Agent(satn).xHatUKF(1,1:i),'r-',...%-Agent(satn).iner_ECI(1,1)*ones(1,length(ttime)),'r-',...
        ttime,Agent(satn).xHatUKF(2,1:i),'b-',...
        ttime, Agent(satn).xHatUKF(3,1:i),'c-','LineWidth',2);
    legend('x','y','z')
    hold on
    %plot(ttime, myEulerCopy(6*(satn-1)+1:6*(satn-1)+3,1:i)','--','LineWidth',2);
    plot(ttime, Agent(satn).iner_ECI(1,1:i),'r--',...%-Agent(satn).iner_ECI(1,1)*ones(1,length(ttime)),'r--',...
        ttime,Agent(satn).iner_ECI(2,1:i)','b--',...
        ttime,Agent(satn).iner_ECI(3,1:i)','c--','LineWidth',2);
    hold off
    grid on
    ylabel('[km]')
    
    ax(2) = subplot(2,1,2);
    plot(ttime,Agent(satn).xHatUKF(4,1:i),'r-',...
        ttime,Agent(satn).xHatUKF(5,1:i),'b-',...
        ttime,Agent(satn).xHatUKF(6,1:i),'c-','LineWidth',2);
    legend('$\dot x$','$\dot y$','$\dot z$')
    hold on
    
    %plot(ttime, myEulerCopy(6*(satn-1)+4:6*(satn-1)+6,1:i)','--','LineWidth',2);
    plot(ttime, Agent(satn).iner_ECI(4,1:i)','r--',...
        ttime,  Agent(satn).iner_ECI(5,1:i)','b--',...
        ttime,  Agent(satn).iner_ECI(6,1:i)','c--','LineWidth',2);
    hold off
    grid on
    xlim([time(1) time(i)])
    ylabel('[km/s]')
    
   
%     % decrease distance between subplots
%     for h = 1:length(ax)
%         vCurrPos = get(ax(h), 'position'); % current position
%         %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
%         set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
%         if h < length(ax)
%             set(ax(h), 'XTickLabel', ' ')
%         end
%     end
     linkaxes(ax,'x');
%     
    
    %%
   figure('Name', ['Agent '  num2str(satn) ': estimation errors'])
   
    ax(1) = subplot(2,1,1);
    plot(ttime,Agent(satn).iner_ECI_EstimationError(1,1:i),'r-',...
        ttime,Agent(satn).iner_ECI_EstimationError(2,1:i),'b-',...
        ttime,Agent(satn).iner_ECI_EstimationError(3,1:i),'c-','LineWidth',2);
    grid on
    ylabel('[m]')
    legend('$E_x$','$E_y$','$E_z$')

    ax(2) = subplot(2,1,2);
    plot(ttime,Agent(satn).iner_ECI_EstimationError(4,1:i),'r-',...
        ttime,Agent(satn).iner_ECI_EstimationError(5,1:i),'b-',...
        ttime,Agent(satn).iner_ECI_EstimationError(6,1:i),'c-','LineWidth',2);
    grid on
    xlim([time(1) time(i)])
    legend('$E_{\dot x}$','$E_{\dot y}$','$E_{\dot z}$')
    ylabel('[m/s]')
  
    % decrease distance between subplots
%     for h = 1:length(ax)
%         vCurrPos = get(ax(h), 'position'); % current position
%         %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
%         set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
%         if h < length(ax)
%             set(ax(h), 'XTickLabel', ' ')
%         end
%     end
     linkaxes(ax,'x');
    
    
    
    
    figure('Name', 'Agents estimation errors')
   
    temp = zeros(ObserverTest.Nagents,i);
    for satn=1:ObserverTest.Nagents,
        temp(satn,:) = Agent(satn).NormEstimationErrorXYZ(1:i);
    end
    
    clear ax
    ax(1) = subplot(3,1,1);
    semilogy(ttime,temp,'LineWidth',2);
    hold on
    semilogy(ttime, norm(ObserverTest.GPSGaussianCovariance(1:3))*ones(1,i),'k--','LineWidth',2);
    hold off
    grid on
    ylabel('$||E_i(x,y,z)||$')
    
    ax(2) = subplot(3,1,2);
    temp = zeros(ObserverTest.Nagents,i);
    for satn=1:ObserverTest.Nagents,
        temp(satn,:) = sqrt(Agent(satn).iner_ECI_EstimationError(4,1:i).^2+Agent(satn).iner_ECI_EstimationError(5,1:i).^2+...
        Agent(satn).iner_ECI_EstimationError(6,1:i).^2);
        %temp(satn,1:i) = Agent(satn).UWBcorrectionActive(1:i)';
    end
    semilogy(ttime,temp,'LineWidth',2)
    grid on
    ylabel('$||dE_i(dx,dy,dz)||$')
    
    
    ax(3) = subplot(3,1,3);
    temp = zeros(ObserverTest.Nagents,i);
    for satn=1:ObserverTest.Nagents,
        temp(satn,1:i) = Agent(satn).MismatchAprioriDistances(1:i)';
    end
    plot(ttime,mean(temp),'-b','LineWidth',2)
    hold on

    temp = zeros(ObserverTest.Nagents,i);
    for satn=1:ObserverTest.Nagents,
        temp(satn,1:i) = Agent(satn).MismatchAposterioriDistances(1:i)';
    end
    plot(ttime,mean(temp),'-r','LineWidth',3)
    
    temp = zeros(ObserverTest.Nagents,i);
    for satn=1:ObserverTest.Nagents,
        temp(satn,1:i) = Agent(satn).ErrorAposterioriDistances(1:i)';
    end
    plot(ttime,mean(temp),':r','LineWidth',3)
    hold off
    grid on
    %ylim([0 10])
    ylabel('SumDistances')
    legend('Pre mismatch','Post mismatch','Post error');
    xlabel('Time (s)')
    
    
%     % decrease distance between subplots
%     for h = 1:length(ax)
%         vCurrPos = get(ax(h), 'position'); % current position
%         %set(ha(h), 'position', (vCurrPos.*[1 1 nF nF])-[vCurrPos(3)*(nF-1)/2 vCurrPos(4)*(nF-1)/2 0 0]);
%         set(ax(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
%         if h < length(ax)
%             set(ax(h), 'XTickLabel', ' ')
%         end
%     end
%     linkaxes(ax,'x');
%     
    

    figure('Name', 'Agents position estimation errors (norm)')
   
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.Mean(ObserverTest.Ntest,:)*1000,ObserverTest.Sigma(ObserverTest.Ntest,:)*1000,'LineWidth',2);
    ylabel('Position: mean+-sigma [m]')
    xlabel('Agents')
    grid on
    
    figure('Name', 'Agents relative distance estimation errors (norm)')
   
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.MeanDistances(ObserverTest.Ntest,:)*1000,ObserverTest.SigmaDistances(ObserverTest.Ntest,:)*1000,'b--','LineWidth',2);
    ylabel('Relative Distances: mean+-sigma [m]')
    xlabel('Agents')
    grid on
    
    if(ObserverTest.UWBoptimizationOn)
        %%
        figure('Name', 'GPS errors with UWB measures')
        
        temp1 = zeros(ObserverTest.Nagents,i);
        temp2 = zeros(ObserverTest.Nagents,i);
        temp3 = zeros(ObserverTest.Nagents,i);
        temp4 = zeros(ObserverTest.Nagents,i);
        temp5 = zeros(ObserverTest.Nagents,i);
        for satn=1:ObserverTest.Nagents,
            temp = Agent(satn).iner_ECI([1 2 3],1:i) - Agent(satn).GPS(:,1:i);
            temp1(satn,:) = sqrt(sum(temp.*temp));
            temp = Agent(satn).iner_ECI([1 2 3],1:i) - Agent(satn).GPSopt(:,1:i);
            temp2(satn,:) = sqrt(sum(temp.*temp));
            temp = Agent(satn).iner_ECI([4 5 6],1:i) - Agent(satn).GPSpeed(:,1:i);
            temp3(satn,:) = sqrt(sum(temp.*temp));
            temp = Agent(satn).iner_ECI([4 5 6],1:i) - Agent(satn).GPSpeedOpt(:,1:i);
            temp4(satn,:) = sqrt(sum(temp.*temp));
            temp = Agent(satn).GPS(:,1:i) - Agent(satn).GPSopt(:,1:i);
            temp5(satn,:) = sqrt(sum(temp.*temp));
        end
        clear ax
        ax(1) = subplot(5,1,1);
        plot(ttime,temp1,'LineWidth',2);
        hold off
        grid on
        ylabel('$E_{GPS}$')
        
        ax(2) = subplot(5,1,2);
        plot(ttime,temp2,'LineWidth',2);
        hold off
        grid on
        ylabel('$E_{GPS,opt}$ ')
        
        ax(3) = subplot(5,1,3);
        semilogy(ttime, sum(temp1),ttime, sum(temp2),'LineWidth',2);
        hold off
        grid on
        legend('$\Sigma ||E_{i,GPS}||$','$\Sigma|| E_{i,GPSopt}||$')
        
        ax(4) = subplot(5,1,4);
        plot(ttime, temp5,'LineWidth',2);
        hold off
        grid on
        legend('$GPS-GPS_{opt}$')
        
        
        ax(5) = subplot(5,1,5);
        semilogy(ttime, sum(temp3),ttime, sum(temp4),'LineWidth',2);
        hold off
        grid on
        legend('$\Sigma ||E_{i,GPSpeed}||$','$\Sigma||E_{i,GPSpeedopt}||$')
        
        xlabel('Time (s)')
        
        
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
        %             %         for satn=1:ObserverTest.Nagents,
        %             %             coldraw = [colore{satn} 'x'];
        %             %             temp = Agent(satn).iner_ECI([1 2 3],kt);
        %             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
        %             %             hold on
        %             %             coldraw = [colore{satn} 'o'];
        %             %             temp = Agent(satn).xHatUKF([1 2 3],kt);
        %             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
        %             %             coldraw = [colore{satn} 's'];
        %             %             temp = Agent(satn).GPS(:,kt);
        %             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
        %             %             coldraw = [colore{satn} 'd'];
        %             %             temp = Agent(satn).GPSopt(:,kt);
        %             %             plot3(temp(1),temp(2),temp(3),coldraw,'MarkerSize',8);
        %             %         end
        %             %         hold off
        %             %         grid on
        %             %         pause
        %             for satn=1:ObserverTest.Nagents,
        %                 ax(1) = subplot(2,1,1);
        %                 coldraw = [colore{satn} 'x'];
        %                 temp = Agent(satn).iner_ECI([1 2 3],kt);
        %                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
        %                 hold on
        %                 temp = Agent(satn).iner_ECI([1 2 3],kt-1:kt);
        %                 plot(temp(1,:),temp(2,:),[colore{satn} '-'],'LineWidth',2);
        %
        %                 coldraw = [colore{satn} 'o'];
        %                 temp = Agent(satn).xHatUKF([1 2 3],kt);
        %                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
        %                 temp = Agent(satn).xHatUKF([1 2 3],kt-1:kt);
        %                 plot(temp(1,:),temp(2,:),[colore{satn} '--'],'LineWidth',2);
        %
        %                 coldraw = [colore{satn} 's'];
        %                 temp = Agent(satn).GPS(:,kt);
        %                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
        %
        %                 coldraw = [colore{satn} 'd'];
        %                 temp = Agent(satn).GPSopt(:,kt);
        %                 plot(temp(1),temp(2),coldraw,'MarkerSize',8);
        %
        %
        %                 ax(2) = subplot(2,1,2);
        %                 coldraw = [colore{satn} 'x'];
        %                 temp = Agent(satn).iner_ECI([1 2 3],kt);
        %                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
        %                 hold on
        %                 temp = Agent(satn).iner_ECI([1 2 3],kt-1:kt);
        %                 plot(temp(1,:),temp(3,:),[colore{satn} '-'],'LineWidth',2);
        %
        %                 coldraw = [colore{satn} 'o'];
        %                 temp = Agent(satn).xHatUKF([1 2 3],kt);
        %                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
        %                 temp = Agent(satn).xHatUKF([1 2 3],kt-1:kt);
        %                 plot(temp(1,:),temp(3,:),[colore{satn} '--'],'LineWidth',2);
        %
        %
        %                 coldraw = [colore{satn} 's'];
        %                 temp = Agent(satn).GPS(:,kt);
        %                 plot(temp(1),temp(3),coldraw,'MarkerSize',8);
        %                 coldraw = [colore{satn} 'd'];
        %                 temp = Agent(satn).GPSopt(:,kt);
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
    
end