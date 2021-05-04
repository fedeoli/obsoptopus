%% plot section for cone analysis
function plot_cone_v1(DynOpt)
    % Attitude
    % EULER ANGLES FORM
    if 1
        offset = DynOpt.integration_pos*6;
        offset_plot = 0;
        figure('Name','Attitude estimate - Euler angles')
        for k=1:3
            n_subplot = 3;
            ax_index = k-offset_plot;
            ax(ax_index)=subplot(n_subplot,1,ax_index);
            hold on
            % results
            plot(DynOpt.time,DynOpt.cone_quat_low(k,:),'-','LineWidth' ,2)
            plot(DynOpt.time,DynOpt.cone_quat_high(k,:),'-','LineWidth' ,2);
            plot(DynOpt.time,DynOpt.Opt_quat_runtime(k,:),'--','LineWidth' ,2);
            % reference
            plot(DynOpt.time(1:length(DynOpt.desatt_true(k,:))),DynOpt.desatt_true(k,:),'r:','LineWidth' ,2);
            ylabel(strcat('X_',num2str(k),' [rad]'));
            xlabel('simulation time [s]')
            grid on
        end
    end

    if 1
        % Rotational velocity
        figure('Name','Rotational velocity estimate')
        offset = DynOpt.integration_pos*6;
        for k=1:3
            n_subplot = 3;
            ax_index = k;
            ax(ax_index)=subplot(n_subplot,1,ax_index);
            hold on
            plot(DynOpt.time,DynOpt.coneState_low(offset+4+k,:),'-','LineWidth' ,2)
            plot(DynOpt.time,DynOpt.coneState_high(offset+4+k,:),'-','LineWidth' ,2);
            plot(DynOpt.time,DynOpt.OptXstory_runtime(offset+4+k,:),'--','LineWidth' ,2);
            ylabel(strcat('w_',num2str(k),' [rad/s]'));
            xlabel('simulation time [s]');
            grid on
        end
    end
end