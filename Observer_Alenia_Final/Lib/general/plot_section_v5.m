%% plot section 
function plot_section_v5(DynOpt,params)

    %% plant test

    if DynOpt.ObserverOn == 0
        if  0 && (DynOpt.integration_pos == 1)
            figure(1)
            hold on
            grid on

            title('Fleet trajectories')
            xlabel('x Km')
            ylabel('y Km')
            zlabel('z Km')

            nagent = params.Nagents;

            % real trajectories
            for i = 1:nagent
                color = [rand rand rand];
                style_real = '.';
                linewidth = 1;

                Chi = DynOpt.position_state(1+6*(i-1):3+6*(i-1),:);
                xreal = Chi(1,:);
                yreal = Chi(2,:);
                zreal = Chi(3,:);
                plot3(xreal,yreal,zreal,style_real,'MarkerFaceColor',color);
                plot3(xreal(1),yreal(1),zreal(1),'bo')
                plot3(xreal(end),yreal(end),zreal(end),'bo')
            end
        end
        
        % Attitude
        % EULER ANGLES FORM
        if 1
            start = 1;
            stop = 3;
            offset = 0;
            figure('Name','Attitude estimate - Euler angles')
            for k=start:stop
                n_subplot = 3;
                ax_index = k-offset;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                % results
                plot(DynOpt.time,DynOpt.True_quat(k,:),':','LineWidth' ,2);
                % reference
                plot(DynOpt.time(1:length(DynOpt.desatt_true(k,:))),DynOpt.desatt_true(k,:),'r--','LineWidth' ,2);
                ylabel(strcat('X_',num2str(k),' [rad]'));
                xlabel('simulation time [s]')
                grid on
                legend('state','ref')
            end
            linkaxes(ax,'x');
        end
    end

    %% estimation 
    if DynOpt.ObserverOn
        % ESTIMATED STATES
        if 1
            % Position
            if 0 && DynOpt.integration_pos == 1
                figure('Name','Position estimate')
                for k=1:3
                    ax_index = k;
                    n_subplot = 3;
                    ax(ax_index)=subplot(n_subplot,1,ax_index);
                    hold on
                    
                    try
                        plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    catch
                        plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    end
                    
                    ylabel(strcat('X_',num2str(k),' [Km]'));
                    xlabel('simulation time [s]')
                    grid on
                    legend('Opt','True')
                end
                linkaxes(ax,'x');

                % Velocity 
                figure('Name','Velocity estimate')
                for k=4:6
                    n_subplot = 3;
                    ax_index = k-3;
                    ax(ax_index)=subplot(n_subplot,1,ax_index);
                    hold on
                    try
                        plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    catch
                        plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    end
                    ylabel(strcat('X_',num2str(k),' [Km]'));
                    xlabel('simulation time [s]')
                    grid on
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end

            % Attitude
            if DynOpt.integration_att == 1

                % QUATERNION FORM
                if DynOpt.integration_pos == 1
                    start = 7;
                    stop = 10;
                    offset = 6;
                else
                    start = 1;
                    stop = 4;
                    offset = 0;
                end
                figure('Name','Attitude estimate - quaternions')
                for k=start:stop
                    n_subplot = 4;
                    ax_index = k-offset;
                    ax(ax_index)=subplot(n_subplot,1,ax_index);
                    hold on
                    % results
                    try
                        plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    catch
                        plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    end
                    ylabel(strcat('X_',num2str(k)));
                    xlabel('simulation time [s]')
                    grid on
                    legend('Opt','True')
                end
                linkaxes(ax,'x');

                % Attitude
                % EULER ANGLES FORM
                if 1
                    start = 1;
                    stop = 3;
                    offset = 0;
                    figure('Name','Attitude estimate - Euler angles')
                    for k=start:stop
                        n_subplot = 3;
                        ax_index = k-offset;
                        ax(ax_index)=subplot(n_subplot,1,ax_index);
                        hold on
                        % results
                        try
                            plot(DynOpt.time,DynOpt.Opt_quat_runtime(k,:)*180/pi,'-',DynOpt.time,DynOpt.True_quat(k,:)*180/pi,':','LineWidth' ,2);
                        catch
                            plot(DynOpt.time,DynOpt.Opt_quat(k,:)*180/pi,'-',DynOpt.time,DynOpt.True_quat(k,:)*180/pi,':','LineWidth' ,2);
                        end
                        % reference
                        plot(DynOpt.time(1:length(DynOpt.desatt_true(k,:))),DynOpt.desatt_true(k,:)*180/pi,'r--','LineWidth' ,2);
                        ylabel(strcat('X_',num2str(k),' [rad]'));
                        xlabel('simulation time [s]')
                        grid on
                        legend('Opt','True')
                    end
                    linkaxes(ax,'x');
                end

                % Rotational velocity
                figure('Name','Rotational velocity estimate')
                if DynOpt.integration_pos == 1
                    start = 11;
                    stop = 13;
                    offset = 10;
                else
                    start = 5;
                    stop = 7;
                    offset = 4;
                end
                for k=start:stop
                    n_subplot = 3;
                    ax_index = k-offset;
                    ax(ax_index)=subplot(n_subplot,1,ax_index);
                    hold on
                    try
                        plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    catch
                        plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    end
                    ylabel(strcat('w_',num2str(k),' [rad/s]'));
                    xlabel('simulation time [s]');
                    grid on
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end

            % Param estimate
            %%% bias
            figure('Name','Param estimate')
            n_param = length(params.param_estimate);
            for k=1:n_param-3
                n_subplot = n_param-3;
                ax_index = k;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                try
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(DynOpt.StateDim+k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((DynOpt.StateDim+k),:),':','LineWidth' ,2);
                catch
                    plot(DynOpt.time,DynOpt.OptXstory(DynOpt.StateDim+k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((DynOpt.StateDim+k),:),':','LineWidth' ,2);
                end
                ylabel(strcat('P_',num2str(k),' [rad/s]'));
                xlabel('simulation time [s]');
                grid on
                legend('Opt','True')
            end
            linkaxes(ax,'x');
            %%% inertias
            figure('Name','Param estimate - I')
            n_param = length(params.param_estimate);
            for k=n_param-2:n_param
                n_subplot = 3;
                ax_index = k-n_param+3;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                try
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(DynOpt.StateDim+k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((DynOpt.StateDim+k),:),':','LineWidth' ,2);
                catch
                    plot(DynOpt.time,DynOpt.OptXstory(DynOpt.StateDim+k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((DynOpt.StateDim+k),:),':','LineWidth' ,2);
                end
                ylabel(strcat('P_',num2str(k),' [Kg m^2]'));
                xlabel('simulation time [s]');
                grid on
                legend('Opt','True')
            end
            linkaxes(ax,'x');
        end


        % WINDOWED DATA
        if 1
            figure('Name','Windowed data')
            for k=1:length(params.observed_state)
                n_subplot = length(params.observed_state);
                ax_index = k;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                WindowTime = DynOpt.time(DynOpt.temp_time);
                OptTime = DynOpt.time(DynOpt.opt_chosen_time);
                plot(WindowTime,DynOpt.Y_full_story(k,DynOpt.temp_time),'s','LineWidth',2,'MarkerSize',5);
                plot(OptTime,DynOpt.Y_full_story(k,DynOpt.opt_chosen_time),'r+','LineWidth',2,'MarkerSize',5);
                try
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(params.observed_state(k),:),'--','LineWidth',2)
                catch
                    plot(DynOpt.time,DynOpt.OptXstory(params.observed_state(k),:),'--','LineWidth',2)
                end
                plot(DynOpt.time,DynOpt.OptXstoryTRUE(params.observed_state(k),:).^DynOpt.measure_exp ,'b-','LineWidth',2)
                ylabel(strcat('w_',num2str(k),' [rad/s]'));
                xlabel('simulation time [s]');
                legend('measures','opt shots','estimation','true')
            end
            
            figure('Name','Windowed data - Magnetometers')
            nMagneto_out = 3*DynOpt.nMagneto;
            for k=1:nMagneto_out
                index = 3+k;
                WindowTime = DynOpt.time(DynOpt.temp_time);
                OptTime = DynOpt.time(DynOpt.opt_chosen_time);
                n_subplot = nMagneto_out;
                ax_index = k;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                plot(WindowTime,DynOpt.Y_full_story(index,DynOpt.temp_time),'s','LineWidth',2,'MarkerSize',5);
                plot(OptTime,DynOpt.Y_full_story(index,DynOpt.opt_chosen_time),'r+','LineWidth',2,'MarkerSize',5);
                plot(DynOpt.Yhat_full_story(index,:),'--','LineWidth',2,'MarkerSize',5);
                plot(DynOpt.Y_full_story(index,:),'LineWidth',2,'MarkerSize',5);
                ylabel(strcat('M',num2str(k)));
                xlabel('simulation time [s]');
                legend('measures','opt shots','estimation','true')
            end
        end
        
        
        % COST FUNCTION
        if 1
            figure('Name','J index')
            J_time = DynOpt.opt_chosen_time;
            semilogy(J_time,DynOpt.Jstory,'b','LineWidth',2)
            xlabel('simulation time [s]')
            ylabel('cost function')
        end
        
        % COMPARE COST FUNCTIONS
        if 0
            window_interval = 5;
            J_time = DynOpt.opt_chosen_time(window_interval:end);
            figure('Name','Cost function comparison')
            hold on 
            grid on
            fplot = @plot;
            fplot(J_time,DynOpt.J_meas(window_interval:end),'LineWidth',1.5)
            fplot(J_time,DynOpt.J_der(window_interval:end),'LineWidth',1.5)
            fplot(J_time,DynOpt.J_int(window_interval:end),'LineWidth',1.5)
            fplot(J_time,DynOpt.J_dyn(window_interval:end),'LineWidth',1.5)
            fplot(J_time,DynOpt.J_quat(window_interval:end),'LineWidth',1.5)
            legend('Meas','Der','Int','Dyn','Quat')
        end

        % DJ THRESHOLDS
        if 1
            figure('Name','dJ condition')
            hold on 
            grid on
            nDj = length(DynOpt.dJ_cond_story(4,:));
            WindowTime = DynOpt.time(DynOpt.temp_time);
            plot(DynOpt.time(1:nDj),DynOpt.dJ_cond_story(4,:));
            plot(DynOpt.time(1:nDj),ones(1,nDj)*DynOpt.dJ_1,'--','LineWidth',2);
            plot(DynOpt.time(1:nDj),ones(1,nDj)*DynOpt.dJ_2,'--','LineWidth',2);
            plot(WindowTime(1:end-1),DynOpt.dJ_cond_story(end,DynOpt.temp_time(1:end-1)),'s','LineWidth',2,'MarkerSize',5);
            legend('dJ condition','dJ thresh')
        end

        % OPTIMISATION TIME
        if 1
            figure('Name','Optimisation times')
            histogram(DynOpt.opt_time);
            xlabel('optimisation time [s]')
            ylabel('occurrences')
        end
        
        % Attitude
        % EULER ANGLES FORM
        if 1
            start = 1;
            stop = 3;
            offset = 0;
            figure('Name','Estimation error - Euler angles ')
            for k=start:stop
                n_subplot = 3;
                ax_index = k-offset;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                % results
                plot(DynOpt.time,DynOpt.OptErrorStory_Euler(k,:)*180/pi,'-','LineWidth' ,2);
                ylabel(strcat('X_',num2str(k),' [deg]'));
                xlabel('simulation time [s]')
                xlim([0;DynOpt.Tend]);
                grid on
            end
%             linkaxes(ax,'x');
        end
    end
end
