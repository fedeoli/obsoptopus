%% Montecarlo simulations 
function struct = grid_analyse(window_interval,window_step,sample_interval,sample_step,pathname,plot_flag)
    
    %% Data analysis
    path = pathname;
    n_window = length(window_interval(1):window_step:window_interval(2));
    n_samples = length(sample_interval(1):sample_step:sample_interval(2));
        
    for counter_window = window_interval(1):window_step:window_interval(2)
        for counter_sample = sample_interval(1):sample_step:sample_interval(2)
            
            file = strcat('/simulation_',int2str(counter_window));
            file = strcat(file,'-');
            file = strcat(file,int2str(counter_sample));
            
            final_path = strcat(path,file);
            load(final_path);
            
            %%%%%% corrections on datasets %%%%%%
            % error definition
%             DynOpt.OptErrorStory = DynOpt.OptXstoryTRUE - DynOpt.OptXstory;
%             save(final_path)
            
            % performance index computation (troppo lungo-rilancia direttamente la simulazione)
%             [DynOpt.lambda_min,~] = performance_index_v3(DynOpt.Y_full_story);
%             save(final_path)
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            row = counter_window/window_step;
            column = counter_sample/sample_step;
            lambda_min(row,column) = DynOpt.lambda_min;
            mean_error(row,column) = norm(sum(abs(DynOpt.OptErrorStory),2));
            
            keep path lambda_min mean_error counter_window window_interval counter_sample sample_interval sample_step window_step n_window n_samples plot_flag  
            temp = strcat(path,'/recap');
            save(temp);            
        end
    end
    
    load(temp)

    struct.lambda_min = real(lambda_min(:,2:end));
    struct.mean_error = real(mean_error(:,2:end));
    
    struct.window_interval = window_interval;
    struct.sample_interval = sample_interval;
    
    [struct.row_lambda, struct.col_lambda] = find(ismember(struct.lambda_min, min(struct.lambda_min(:))));
    [struct.row_error, struct.col_error] = find(ismember(struct.mean_error, min(struct.mean_error(:))));
    
    if plot_flag == 1
        figure
        title('Lambda (MIN)')
        hold on
        surf(struct.lambda_min)
        h1 = plot3(struct.col_lambda,struct.row_lambda,struct.lambda_min(struct.row_lambda,struct.col_lambda),'ro');
        set(h1, 'markerfacecolor', get(h1, 'color'));
        
        figure
        title('ERROR (MIN)')
        hold on
        surf(struct.mean_error)
        h1 = plot3(struct.col_error,struct.row_error,struct.mean_error(struct.row_error,struct.col_error),'ro');       
        set(h1, 'markerfacecolor', get(h1, 'color'));
    end
        
end