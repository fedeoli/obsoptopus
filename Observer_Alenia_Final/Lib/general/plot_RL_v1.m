%% plot section 
function plot_RL_v1(RL,s)

    % get number of iterations - RL
    setup = RL.S.setup;
    Nplot = length(RL.S.search{s}.DynOpt);
    time_total = [];
    
    % plot flags
    flag_vel_att = 1;
    flag_att = 1;
    flag_eul = 0;
    
    %%% define here all the figures
    %%%%% angular velocity %%%%%
    if flag_vel_att
        vel = figure('Name','Velocity estimate');
        for k=1:3
            n_subplot = 3;
            ax(k)=subplot(n_subplot,1,k);
            ylabel(strcat('X_',num2str(k),' [rad/s]'));
            xlabel('simulation time [s]')
            grid on
        end
        linkaxes(ax,'x');
    end
    
    %%%%% quaternions %%%%%
    if flag_att
        quat = figure('Name','Quaternion estimate');
        for k=1:4
            n_subplot = 4;
            ax(k)=subplot(n_subplot,1,k);
            ylabel(strcat('q_',num2str(k)));
            xlabel('simulation time [s]')
            grid on
        end
        linkaxes(ax,'x');
    end
    
    %%%%% quaternions %%%%%
    if flag_eul
        eul = figure('Name','Euler angles estimate');
        for k=1:3
            n_subplot = 3;
            ax(k)=subplot(n_subplot,1,k);
            ylabel(strcat('x_',num2str(k)));
            xlabel('simulation time [s]')
            grid on
        end
        linkaxes(ax,'x');
    end
    
    for i=1:Nplot
        % assign DynOpt
        DynOpt = RL.S.search{s}.DynOpt{i};
        if i>1
            time = DynOpt.time;
            time_total = [time_total, time];
        else
            time = DynOpt.time;
            time_total = time;
        end
        
        if setup.ObserverOn
            % ESTIMATED STATES
            
            % Angular velocity 
            if flag_vel_att
                figure(vel);
                handles=flipud(findobj(vel,'Type','axes'));
                xlim([time_total(1) time_total(end)]);
                hold on
                for k=1:3
                    offset_pos = setup.integration_pos*6;
                    offset = offset_pos + 4;
                    axes(handles(k));
                    hold on
                    plot(time,DynOpt.OptXstory_runtime(offset+k,:),'r-',time,DynOpt.OptXstoryTRUE((offset+k),:),'b:','LineWidth' ,2);
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end
            
            % Quaternions
            if flag_att
                figure(quat);
                handles=flipud(findobj(quat,'Type','axes'));
                xlim([time_total(1) time_total(end)]);
                hold on
                for k=1:4
                    offset_pos = setup.integration_pos*6;
                    offset = offset_pos;
                    axes(handles(k));
                    hold on
                    plot(time,DynOpt.OptXstory_runtime(offset+k,:),'r-',time,DynOpt.OptXstoryTRUE((offset+k),:),'b:','LineWidth' ,2);
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end
            
            % Euler angles 
            if flag_eul
                figure(eul);
                handles=flipud(findobj(eul,'Type','axes'));
                xlim([time_total(1) time_total(end)]);
                hold on
                for k=1:3
                    offset = 0;
                    axes(handles(k));
                    hold on
                    plot(time,DynOpt.Opt_quat_runtime(offset+k,:),'r-',time,DynOpt.True_quat((offset+k),:),'b:','LineWidth' ,2);
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end
        end
    end

end
