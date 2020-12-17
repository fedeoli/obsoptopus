%% plot section 

global DynOpt params

%% plant test
if 1
    figure()
    if DynOpt.simulationModel == 1
        for i=1:DynOpt.StateDim
            ax(i)=subplot(DynOpt.StateDim,1,i);
            plot(DynOpt.time,DynOpt.stateStory(i,:),'LineWidth',2);
            grid on
            ylabel(strcat('x_',int2str(i)));
            title('Simulation test ')
            legend(strcat('x_',int2str(i)))
        end
    else
        plot(DynOpt.time,DynOpt.stateStory(1,:),'LineWidth',2);
        grid on
        ylabel(strcat('x_',int2str(1)));
        title('Simulation test ')
        legend(strcat('x_',int2str(1)))
    end

    if strcmp(DynOpt.modelname,'runaway') & (DynOpt.simulationModel == 1)
       figure()
       hold on
       grid on
       plot(DynOpt.stateStory(1,:),DynOpt.stateStory(2,:),'b--') 
       plot(DynOpt.stateStory(1,1),DynOpt.stateStory(2,1),'ro')
       plot(DynOpt.stateStory(1,end),DynOpt.stateStory(2,end),'r+')
       legend('evolution','init','end')
    end
end

%% estimation 
if DynOpt.ObserverOn
    if 1
    figure('Name','J index')
    J_time = linspace(0,DynOpt.Tend,length(DynOpt.Jstory));
    semilogy(J_time,DynOpt.Jstory,'b','LineWidth',2)
    end

    %%% CON IL WRONG %%%
    if 1
        if DynOpt.simulationModel == 1
            figure('Name','Estimates shots')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                hold on
                plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                ylabel(strcat('X_',num2str(k)));
                plot(DynOpt.time,DynOpt.Xstory(k,:),'r--');
                grid on
                legend('Opt','True','Wrong')
            end
            linkaxes(ax,'x');
        else
            figure('Name','Estimates shots')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                hold on
                plot(DynOpt.time,DynOpt.OptXstory(k,:),'-','LineWidth' ,2);
                if k == DynOpt.params.observed_state
                    plot(DynOpt.time,DynOpt.stateStory,':','LineWidth' ,2)
                    ylabel(strcat('X_',num2str(k)));
                    grid on
                    legend('Opt','True')
                else
                    grid on
                    legend('Opt')
                end
            end
            linkaxes(ax,'x');
        end
    end
    
    %%% SENZA IL WRONG %%%
    if 0
        if DynOpt.simulationModel == 1
            figure('Name','Estimates shots ')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                hold on
                plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                ylabel(strcat('X_',num2str(k)));
                grid on
                legend('Opt','True')
            end
            linkaxes(ax,'x');
        else
            figure('Name','Estimates shots ')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                hold on
                plot(DynOpt.time,DynOpt.OptXstory(k,:),'-','LineWidth' ,2);
                if k == DynOpt.params.observed_state
                    plot(DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2)
                    grid on
                    legend('Opt','True')
                else
                    grid on
                    legend('Opt')
                end
            end
            linkaxes(ax,'x');
        end
    end

    if 0
        if DynOpt.simulationModel == 1
            figure('Name','Estimation Errors')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                plot(DynOpt.time,abs(DynOpt.OptErrorStory(k,:)),'LineWidth',2);
                ylabel(strcat('E_',num2str(k)));
                grid on
            end
            linkaxes(ax,'x');
        end
    end

    %% windowed data - con il wrong
    if 0
        figure('Name','Check windowed data')
        clear ax
        hold on
        WindowTime = DynOpt.time(DynOpt.temp_time);
        plot(WindowTime,DynOpt.Y_full_story(1,:),'s','LineWidth',2,'MarkerSize',5);
        plot(DynOpt.time,DynOpt.OptXstory(params.observed_state,:),'--','LineWidth',2)
        if DynOpt.simulationModel == 1
            plot(DynOpt.time,DynOpt.OptXstoryTRUE(params.observed_state,:).^DynOpt.measure_exp ,'b-','LineWidth',2)
            plot(DynOpt.time,DynOpt.Xstory(params.observed_state,:),':','LineWidth',2)
            legend('measures','estimation','true','wrong')
        else
            plot(DynOpt.time,DynOpt.stateStory,':','LineWidth' ,2)
            ylabel(strcat('X_obs'));
            grid on
            legend('measures','estimation','data')
        end
    end
    
    %% windowed data - senza il wrong
    if 1
        figure('Name','Check windowed data')
        clear ax
        hold on
        WindowTime = DynOpt.time(DynOpt.temp_time);
        plot(WindowTime,DynOpt.Y_full_story(1,:),'s','LineWidth',2,'MarkerSize',5);
        plot(DynOpt.time,DynOpt.OptXstory(params.observed_state,:),'--','LineWidth',2)
        if DynOpt.simulationModel == 1
            plot(DynOpt.time,DynOpt.OptXstoryTRUE(params.observed_state,:).^DynOpt.measure_exp ,'b-','LineWidth',2)
            legend('measures','estimation','true')
        else
            plot(DynOpt.time,DynOpt.stateStory,':','LineWidth' ,2)
            ylabel(strcat('X_obs'));
            grid on
            legend('measures','estimation','data')
        end
    end
    
    %% gradient story
    if 0
        figure('Name','Gradient story')
        for k=1:length(DynOpt.X)
            ax(k)=subplot(length(DynOpt.X),1,k);
            plot(DynOpt.grad_story(k,:),'LineWidth',2);
            ylabel(strcat('grad_',num2str(k)));
            grid on
        end
        linkaxes(ax,'x');
    end
    
end
