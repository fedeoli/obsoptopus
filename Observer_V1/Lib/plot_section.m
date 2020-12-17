%% plot section 

global DynOpt

%% plant test
if 1
figure(1)
ax(1)=subplot(4,1,1);
plot(DynOpt.time,DynOpt.stateStory(1,:),'LineWidth',2);
grid on
ylabel('$x_1$','Interpreter','Latex')
title('Simulation test ')
legend('x_1')
if forwardSimulation == 0
    set(gca, 'XDir','reverse')
end

ax(2)=subplot(4,1,2);
plot(DynOpt.time,DynOpt.stateStory(2,:),'LineWidth',2);
grid on
ylabel('$x_2$','Interpreter','Latex')
legend('x_2')
if forwardSimulation == 0
    set(gca, 'XDir','reverse')
end

ax(3)=subplot(4,1,3);
plot(DynOpt.time,DynOpt.stateStory(3,:),'LineWidth',2);
grid on
ylabel('$x_3$','Interpreter','Latex')
legend('x_3')
if forwardSimulation == 0
    set(gca, 'XDir','reverse')
end

ax(4)=subplot(4,1,4);
plot(DynOpt.time,DynOpt.stateStory(4,:),'-','LineWidth',2);
grid on
ylabel('x_4')
linkaxes(ax,'x');
legend('x_4')
if forwardSimulation == 0
    set(gca, 'XDir','reverse')
end
end

%% estimation 
if 1
figure('Name','J index')
semilogy(DynOpt.time,DynOpt.Jstory,'ro')

figure('Name','Estimates shots ')
for k=1:length(DynOpt.X)
    ax(k)=subplot(length(DynOpt.X),1,k);
    eval(['plot(DynOpt.time,DynOpt.OptXstory(' num2str(k) ',:),''-'',DynOpt.time,DynOpt.OptXstoryTRUE(' num2str(k) ',:),'':'',''LineWidth'' ,2)']);
    eval(['ylabel(''X_' num2str(k) ''')' ]);
    grid on
end
linkaxes(ax,'x');

figure('Name','Estimation Errors')
for k=1:length(DynOpt.X)
    ax(k)=subplot(length(DynOpt.X),1,k);
    eval(['plot(DynOpt.time,DynOpt.OptXstoryTRUE(' num2str(k) ',:)-DynOpt.OptXstory(' num2str(k) ',:),''LineWidth'' ,2)']);
    eval(['ylabel(''E_' num2str(k) ''')' ]);
    grid on
end
linkaxes(ax,'x');
end

%% windowed data
if 0
figure('Name','Check windowed data')
clear ax
WindowTime = DynOpt.time(temp_time(end-DynOpt.w+1:1:end));
ax(1) = subplot(2,1,1);
plot(WindowTime,DynOpt.Y(1,:),'s',DynOpt.time,DynOpt.state(2,:),'--','LineWidth',2,'MarkerSize',5)
end
