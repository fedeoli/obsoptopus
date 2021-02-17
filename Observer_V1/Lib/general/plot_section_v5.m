%% plot section 

global DynOpt params

%% plant test
if DynOpt.ObserverOn == 0
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
        plot(DynOpt.time,DynOpt.dataStory(1,:),'LineWidth',2);
        grid on
        ylabel(strcat('x_',int2str(1)));
        title('Simulation test ')
        legend(strcat('x_',int2str(1)))
    end

    if strcmp(DynOpt.modelname,'runaway') && (DynOpt.simulationModel == 1)
       % state trajectory
       figure()
       hold on
       grid on
       plot(DynOpt.stateStory(1,:),DynOpt.stateStory(2,:),'b--') 
       plot(DynOpt.stateStory(1,1),DynOpt.stateStory(2,1),'ro')
       plot(DynOpt.stateStory(1,end),DynOpt.stateStory(2,end),'r+')
       legend('evolution','init','end')
       
       % output
       figure()
       hold on
       grid on
       plot(DynOpt.time,DynOpt.outputStory,'b-','LineWidth',2) 
       if DynOpt.check == 1
           plot(DynOpt.time,DynOpt.dataStory(1,:),'r--','LineWidth',2);
       end
       legend('output')
    end
end

%% estimation 
if DynOpt.ObserverOn
    if 1
    figure('Name','J index')
    J_time = linspace(0,DynOpt.Tend,length(DynOpt.Jstory));
    semilogy(J_time,DynOpt.Jstory,'b','LineWidth',2)
    xlabel('simulation time [s]')
    ylabel('cost function')
    end
    
    %%% SENZA IL WRONG %%%
    if 1
        if DynOpt.simulationModel == 1
            figure('Name','Estimates shots ')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                hold on
                plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                ylabel(strcat('X_',num2str(k)));
                xlabel('simulation time [s]')
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
                    plot(DynOpt.time,DynOpt.OptXstoryTRUE(1,:),':','LineWidth' ,2)
                    grid on
                    legend('Opt','True')
                else
                    grid on
                    legend('Opt')
                end
                xlabel('simulation time [s]')
            end
            linkaxes(ax,'x');
        end
    end

    if 1
        if DynOpt.simulationModel == 1
            figure('Name','Estimation Errors')
            for k=1:length(DynOpt.X)
                ax(k)=subplot(length(DynOpt.X),1,k);
                plot(DynOpt.time,abs(DynOpt.OptErrorStory(k,:)),'LineWidth',2);
                ylabel(strcat('E_',num2str(k)));
                xlabel('simulation time [s]')
                grid on
            end
            linkaxes(ax,'x');
        end
    end

    % WINDOWED DATA
    if 1
        figure('Name','Windowed data')
        for k=1:DynOpt.dim_out
            n_subplot = DynOpt.dim_out;
            ax_index = k;
            ax(ax_index)=subplot(n_subplot,1,ax_index);
            hold on
            WindowTime = DynOpt.time(DynOpt.temp_time);
            OptTime = DynOpt.time(DynOpt.opt_chosen_time);
            plot(WindowTime,DynOpt.Y_full_story(k,DynOpt.temp_time),'s','LineWidth',2,'MarkerSize',5);
            plot(OptTime,DynOpt.Y_full_story(k,DynOpt.opt_chosen_time),'r+','LineWidth',2,'MarkerSize',5);
            plot(DynOpt.time,DynOpt.Y_full_story(k,:),'--','LineWidth',2)
            if DynOpt.simulationModel == 1
                plot(DynOpt.time,DynOpt.Yhat_full_story(k,:).^DynOpt.measure_exp ,'b-','LineWidth',2)
            else
                plot(DynOpt.time,DynOpt.Yhat_full_story(k,:).^DynOpt.measure_exp ,'b-','LineWidth',2)
            end
            ylabel(strcat('X_',num2str(k)));
            xlabel('simulation time [s]')
            legend('measures','chosen','true','estimation')
        end
    end
    
    % CHECK OUTPUT DERIVATIVE
    if 0
        % measures numeric computations
        figure('Name','Measure numeric derivative & integral')
        for k=1:2
            n_subplot = 2;
            ax_index = k;
            ax(ax_index)=subplot(n_subplot,1,ax_index);
            hold on
            plot(DynOpt.time,DynOpt.Y_full_story(k,:),'LineWidth',2);
            plot(DynOpt.time,DynOpt.dY_full_story(k,:),'LineWidth',2);
            plot(DynOpt.time,DynOpt.Y_int(k,:),'o','LineWidth',1);
            legend('Y','dY','Y_int')
        end
        
        % estimations - analytic computations
        figure('Name','Gyro analytic derivative & integral')
        for k=1:2
            n_subplot = 2;
            ax_index = k;
            ax(ax_index)=subplot(n_subplot,1,ax_index);
            hold on
            plot(DynOpt.time,DynOpt.Yhat_full_story(k,:),'LineWidth',2);
            plot(DynOpt.time,DynOpt.dYhat_full_story(k,:),'LineWidth',2);
            plot(DynOpt.time,DynOpt.Yhat_int(k,:),'o','LineWidth',1);
            legend('Y','dY','Y_int')
        end
    end
    
    %% gradient story
    if 1
        figure('Name','Gradient story')
        for k=1:length(DynOpt.X)
            ax(k)=subplot(length(DynOpt.X),1,k);
            plot(DynOpt.grad_story(k,:),'LineWidth',2);
            ylabel(strcat('grad_',num2str(k)));
            xlabel('simulation time [s]')
            grid on
        end
        linkaxes(ax,'x');
    end
    
    %% OPTIMISATION TIME
    if 1
        figure('Name','Optimisation times')
        histogram(DynOpt.opt_time);
    end
end
