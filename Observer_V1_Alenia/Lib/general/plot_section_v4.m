%% plot section 
function plot_section_v4(DynOpt,params,struct)

    %% plant test

    if 1 && (DynOpt.integration_pos == 1)
        figure
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
            linewidth = 2;

            Chi = DynOpt.position_state(1+6*(i-1):3+6*(i-1),:);
            xreal = Chi(1,:);
            yreal = Chi(2,:);
            zreal = Chi(3,:);
            plot3(xreal,yreal,zreal,style_real,'MarkerFaceColor',color);
            plot3(xreal(1),yreal(1),zreal(1),'bo')
            plot3(xreal(end),yreal(end),zreal(end),'bo')
        end
    end

    %% estimation 
    if DynOpt.ObserverOn

        % COST FUNCTION
        if 1
            figure('Name','J index')
            J_time = linspace(0,DynOpt.Tend,length(DynOpt.Jstory));
            semilogy(J_time,DynOpt.Jstory,'b','LineWidth',2)
            xlabel('simulation time [s]')
            ylabel('cost function')
        end

        % ESTIMATED STATES
        if 1
            % Position
            if DynOpt.integration_pos == 1
                figure('Name','Position estimate')
                for k=1:3
                    ax_index = k;
                    n_subplot = 3;
                    ax(ax_index)=subplot(n_subplot,1,ax_index);
                    hold on
                    plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
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
                    plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
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
                    plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
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
                        plot(DynOpt.time,DynOpt.Opt_quat(k,:),'-',DynOpt.time,DynOpt.True_quat(k,:),':','LineWidth' ,2);
                        % reference
    %                     plot(DynOpt.time,DynOpt.desatt_true(k,:),'r--','LineWidth' ,2);
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
                    plot(DynOpt.time,DynOpt.OptXstory(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    ylabel(strcat('w_',num2str(k),' [rad/s]'));
                    xlabel('simulation time [s]');
                    grid on
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end

            % Param estimate
            figure('Name','Param estimate')
            n_param = length(params.param_estimate);
            for k=1:n_param
                n_subplot = n_param;
                ax_index = k;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                plot(DynOpt.time,DynOpt.OptXstory(DynOpt.StateDim+k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((DynOpt.StateDim+k),:),':','LineWidth' ,2);
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
                plot(DynOpt.time,DynOpt.OptXstory(params.observed_state(k),:),'--','LineWidth',2)
                plot(DynOpt.time,DynOpt.OptXstoryTRUE(params.observed_state(k),:).^DynOpt.measure_exp ,'b-','LineWidth',2)
                ylabel(strcat('w_',num2str(k),' [rad/s]'));
                xlabel('simulation time [s]');
                legend('measures','opt shots','estimation','true')
            end
        end

        % CHECK OUTPUT DERIVATIVE
        if 0
            % measures numeric computations
            figure('Name','Gyro numeric derivative & integral')
            for k=1:3
                n_subplot = 3;
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
            for k=1:3
                n_subplot = 3;
                ax_index = k;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                plot(DynOpt.time,DynOpt.Yhat_full_story(k,:),'LineWidth',2);
                plot(DynOpt.time,DynOpt.dYhat_full_story(k,:),'LineWidth',2);
                plot(DynOpt.time,DynOpt.Yhat_int(k,:),'o','LineWidth',1);
                legend('Y','dY','Y_int')
            end
        end

        % OPTIMISATION TIME
        if 1
            figure('Name','Optimisation times')
            histogram(DynOpt.opt_time);
            xlabel('optimisation time [s]')
            ylabel('occurrences')
        end
    end
end
