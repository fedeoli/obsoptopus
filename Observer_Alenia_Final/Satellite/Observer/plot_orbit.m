%% plot orbits - from official data
global ObserverTest

%% GPS OPTIMIZATION 
if 1
figure
hold on
grid on

title('Fleet trajectories')
xlabel('x Km')
ylabel('y Km')
zlabel('z Km')
    
% nagent = ObserverTest.Nagents;
nagent = ObserverTest.Nagents;
% nagent = 1;
    
    % real trajectories
    for i = 1:nagent
        color = [rand rand rand];
        style_real = '.';
        linewidth = 2;
        
        Chi = ObserverTest.AllSimulation(ObserverTest.Ntest).iner_ECI(1+6*(i-1):3+6*(i-1),:);
        xreal = Chi(1,:);
        yreal = Chi(2,:);
        zreal = Chi(3,:);
        plot3(xreal,yreal,zreal,style_real,'MarkerFaceColor',color);
        plot3(xreal(1),yreal(1),zreal(1),'bo')
        plot3(xreal(end),yreal(end),zreal(end),'bo')
    end
    
    % estimated trajectories
    for i = 1:nagent
        color = [rand rand rand];
        style_est = '--';
        linewidth = 2;
        
        Chi_est = ObserverTest.AllSimulation(ObserverTest.Ntest).estimated_iner_ECI(1+6*(i-1):3+6*(i-1),:);
        xest = Chi_est(1,:);
        yest = Chi_est(2,:);
        zest = Chi_est(3,:);
        plot3(xest,yest,zest,style_est,'MarkerFaceColor',color);
        plot3(xest(1),yest(1),zest(1),'+r')
        plot3(xest(end),yest(end),zest(end),'+r')
        
%         xGPS = Agent(i).GPS(1,:);
%         yGPS = Agent(i).GPS(2,:);
%         zGPS = Agent(i).GPS(3,:);
%         plot3(xGPS,yGPS,zGPS,'+r');
    end
end

%% old agents position error
if 0
figure
sgtitle("Agents position estimation error");
nagent = ObserverTest.Nagents;
% nagent = 4;
% nplot = 1;
for i = 1:nagent
    subplot(nagent,1,i)
    grid on;
    hold on
    
    for z = 1:ObserverTest.Npassi
       Chi_norm(z) = norm(Chi(:,z));
       Chi_est_norm(z) = norm(Chi_est(:,z));
       Chi_err(z) = norm(Chi(:,z)-Chi_est(:,z));
    end

    traj_error(:,i) = Chi_err;
    plot(1:ObserverTest.Npassi,traj_error(:,i),'b--','LineWidth',2);
end
end

%% EKF
if 0
% nplot = ObserverTest.Nagents;
nagent = 1;
for agent = 1:nagent
    figure
    lab = ["Roll","Pitch","Yaw"];
    for z=1:3
        subplot(3,1,z)
        grid on;
        hold on
        % True attitude
        plot(1:ObserverTest.Npassi, Agent(agent).q_Euler(:,z)*180/pi, 'b-');
        % True attitude
        plot(1:ObserverTest.Npassi, Agent(agent).q_est_Euler(:,z)*180/pi, 'r--'); 
        legend(lab(z),strcat(lab(z), ' estimated'));
    end
    sgtitle("Agent " + agent + " : RPY attitude estimation");
end

% nplot = ObserverTest.Nagents;
nagent = 1;
for agent = 1:nagent
    figure
    lab = ["Roll","Pitch","Yaw"];
    for z=1:3
        subplot(3,1,z)
        grid on;
        hold on
%         plot((Agent(agent).q_Euler(:,z)-Agent(agent).q_est_Euler(:,z))*180/pi)
        plot((Agent(agent).delq_Euler(:,z))*180/pi)
        legend(strcat(lab(z), ' error'));
    end
    sgtitle("Agent " + agent + " : RPY attitude error");
end
end

if 0
nagent = ObserverTest.Nagents;
for agent = 1:nagent
    figure
    lab = ["X axis","Y axis","Z axis"];
    for z=1:3
        subplot(3,1,z)
        grid on;
        hold on
        xlabel('Time (s)')
        ylabel(lab(z));
        % Estimation error
        plot(ObserverTest.window_interval, Agent(agent).iner_ECI_EstimationError(z,ObserverTest.window_interval), 'b-','LineWidth',2);
        
        % Init GPS error
        temp = reshape(ObserverTest.Gpserror(agent,z,ObserverTest.window_interval),1,length(ObserverTest.window_interval));
        plot(ObserverTest.window_interval, temp, 'r-','LineWidth',2);
        legend('FIlter','GPS')
    end
    sgtitle("Agent " + agent + " : GPS VS estimation error");
end
end
