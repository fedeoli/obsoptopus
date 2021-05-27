%% plot section 
function plot_RL_v1(RL,s,flags)

    % get number of iterations - RL
    setup = RL.S.setup;
    Nplot = length(RL.S.search{s}.DynOpt);
    time_total = [];
    
    % plot flags
    flag_vel_att = flags(1);
    flag_att = flags(2);
    flag_eul = flags(3);
    flag_eul_err = flags(4);
    flag_par = flags(5);
    
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
    
    %%%%% euler angles %%%%%
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

    %%%%% euler angles - error %%%%%
    if flag_eul_err
        eul_err = figure('Name','Euler angles error');
        for k=1:3
            n_subplot = 3;
            ax(k)=subplot(n_subplot,1,k);
            ylabel(strcat('x_',num2str(k)));
            xlabel('simulation time [s]')
            grid on
        end
        linkaxes(ax,'x');
    end
    
    %%%%% params %%%%%
    if flag_par
        par = figure('Name','Params estimation');
        for k=1:3
            n_subplot = 3;
            ax(k)=subplot(n_subplot,1,k);
            ylabel(strcat('P_',num2str(k)));
            xlabel('simulation time [s]')
            grid on
        end
        linkaxes(ax,'x');
    end
    
    i = 1;
    while i<=Nplot
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
                    interval = DynOpt.time(1)+1:DynOpt.time(end)+1;
                    plot(time,DynOpt.OptXstory_runtime(offset+k,interval),'r-',time,DynOpt.OptXstoryTRUE((offset+k),interval),'b:','LineWidth' ,2);
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
                    interval = DynOpt.time(1)+1:DynOpt.time(end)+1;
                    plot(time,DynOpt.OptXstory_runtime(offset+k,interval),'r-',time,DynOpt.OptXstoryTRUE((offset+k),interval),'b:','LineWidth' ,2);
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end
            
            % Euler angles - error 
            if flag_eul_err
                figure(eul_err);
                handles=flipud(findobj(eul_err,'Type','axes'));
                xlim([time_total(1) time_total(end)]);
                hold on
                for k=1:3
                    offset = 0;
                    axes(handles(k));
                    hold on
                    plot(time,DynOpt.OptErrorStory_Euler(offset+k,:)*180/pi,'r-','LineWidth' ,2);
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
                    plot(time,DynOpt.Opt_quat_runtime(offset+k,:)*180/pi,'r-',time,DynOpt.True_quat((offset+k),:)*180/pi,'b:','LineWidth' ,2);
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end
            
            % Params
            if flag_par
                figure(par);
                handles=flipud(findobj(par,'Type','axes'));
                xlim([time_total(1) time_total(end)]);
                hold on
                for k=1:3
                    axes(handles(k));
                    interval = DynOpt.time(1)+1:DynOpt.time(end)+1;
                    hold on
                    plot(time,DynOpt.OptXstory_runtime(end+1-k,interval),'r-',time,DynOpt.OptXstoryTRUE((end+1-k),interval),'b:','LineWidth' ,2);
                    legend('Opt','True')
                end
                linkaxes(ax,'x');
            end
        end
        
        fail_flag = strcmp(RL.S.search{s}.DynOpt{min(i+1,Nplot)},'failed run');
        if fail_flag
           i = Nplot+1; 
        else
           i = i+1;
        end
    end

end
