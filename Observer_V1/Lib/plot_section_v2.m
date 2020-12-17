%% plot section 

global DynOpt params

%% plant test
if 0
figure(1)
ax(1)=subplot(4,1,1);
plot(DynOpt.time,DynOpt.stateStory(1,:),'LineWidth',2);
grid on
ylabel('$x_1$','Interpreter','Latex')
title('Simulation test ')
legend('x_1')

ax(2)=subplot(4,1,2);
plot(DynOpt.time,DynOpt.stateStory(2,:),'LineWidth',2);
grid on
ylabel('$x_2$','Interpreter','Latex')
legend('x_2')

ax(3)=subplot(4,1,3);
plot(DynOpt.time,DynOpt.stateStory(3,:),'LineWidth',2);
grid on
ylabel('$x_3$','Interpreter','Latex')
legend('x_3')

ax(4)=subplot(4,1,4);
plot(DynOpt.time,DynOpt.stateStory(4,:),'-','LineWidth',2);
grid on
ylabel('x_4')
linkaxes(ax,'x');
legend('x_4')
end

%% estimation 
if 0
figure('Name','J index')
semilogy(DynOpt.time,DynOpt.Jstory,'b+','LineWidth',2)
end

if 1
figure('Name','Estimates shots ')
for k=1:length(DynOpt.X)
    ax(k)=subplot(length(DynOpt.X),1,k);
    hold on
    plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
    ylabel(strcat('X_',num2str(k)));
    plot(DynOpt.time,DynOpt.Xstory(k,:),'r--');
    grid on
end
linkaxes(ax,'x');
end

if 0
figure('Name','Estimation Errors')
for k=1:length(DynOpt.X)
    ax(k)=subplot(length(DynOpt.X),1,k);
    plot(DynOpt.time,DynOpt.OptXstoryTRUE(k,:)-DynOpt.OptXstory(k,:),'LineWidth',2);
    ylabel(strcat('E_',num2str(k)));
    grid on
end
linkaxes(ax,'x');
end

%% windowed data
if 1
figure('Name','Check windowed data')
clear ax
hold on
WindowTime = DynOpt.time(DynOpt.temp_time);
plot(WindowTime,DynOpt.Y_story(1,:),'s',DynOpt.time,DynOpt.OptXstory(params.observed_state,:),'--','LineWidth',2,'MarkerSize',5)
plot(DynOpt.time,DynOpt.OptXstoryTRUE(params.observed_state,:),'b-')
legend('measures','estimation','true')
end
