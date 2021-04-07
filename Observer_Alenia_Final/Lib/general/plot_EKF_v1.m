%% plot section 
function plot_EKF_v1(DynOpt,params,struct)

    %% plant test

    if 0 && (DynOpt.integration_pos == 1)
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
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
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
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
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
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
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
                    plot(DynOpt.time,DynOpt.OptXstory_runtime(k,:),'-',DynOpt.time,DynOpt.OptXstoryTRUE((k),:),':','LineWidth' ,2);
                    ylabel(strcat('w_',num2str(k),' [rad/s]'));
                    xlabel('simulation time [s]');
                    grid on
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end

            figure('Name','Magnetometers')
            nMagneto_out = 3*DynOpt.nMagneto;
            for k=1:nMagneto_out
                index = 3+k;
                n_subplot = nMagneto_out;
                ax_index = k;
                ax(ax_index)=subplot(n_subplot,1,ax_index);
                hold on
                plot(DynOpt.time,DynOpt.Yhat_full_story(index,:),'--','LineWidth',2,'MarkerSize',5);
                plot(DynOpt.time,DynOpt.Y_full_story(index,:),'LineWidth',2,'MarkerSize',5);
                ylabel(strcat('M',num2str(k)));
                xlabel('simulation time [s]');
                legend('estimation','true')
            end
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
                grid on
            end
            linkaxes(ax,'x');
        end
    end
end
