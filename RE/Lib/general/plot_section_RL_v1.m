%% plot section 
function plot_section_RL_v1(DynOpt,params,struct)

    %% plant test

    %% estimation 
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
            plot(DynOpt.time,DynOpt.attitude_state((k),:),':','LineWidth' ,2);
            ylabel(strcat('X_',num2str(k)));
            xlabel('simulation time [s]')
            grid on
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
                plot(DynOpt.time,DynOpt.True_quat(k,:),':','LineWidth' ,2);
                ylabel(strcat('X_',num2str(k),' [rad]'));
                xlabel('simulation time [s]')
                grid on
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
            plot(DynOpt.time,DynOpt.attitude_state((k),:),':','LineWidth' ,2);
            ylabel(strcat('w_',num2str(k),' [rad/s]'));
            xlabel('simulation time [s]');
            grid on
        end
        linkaxes(ax,'x');
    end
    
    %%%% OBSERVATIVITY ANALYSIS %%%%
    figure('Name','eigenvalue story')
    plot(DynOpt.time,DynOpt.dtheta_eig_story,'bo')
    
    figure('Name','rank story')
    plot(DynOpt.time,DynOpt.dtheta_rank_story,'bo')
end
