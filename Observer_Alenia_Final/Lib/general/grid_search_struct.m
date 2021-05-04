%% Montecarlo simulations 
function grid_search_struct(window_interval,window_step,sample_interval,sample_step,pathname)
    path = pathname;

    for counter_window = window_interval(1):window_step:window_interval(2)
        for counter_sample = sample_interval(1):sample_step:sample_interval(2)
            file = strcat('/simulation_',int2str(counter_window));
            file = strcat(file,'-');
            file = strcat(file,int2str(counter_sample));
            
            final_path = strcat(path,file);
            save temp

            filename = strcat(final_path,'.mat');
            if ~isfile(filename)
                struct = init_struct_analyse(counter_window,counter_sample);

                [DynOpt, params] = MainOpt_DEZ_general_v6_fun(struct);
                save(final_path);
                clear
                load temp
            else
                delete(filename)
                if counter_sample ~= sample_interval(1)
                   counter_sample = counter_sample - sample_step;
                else
                   counter_sample = sample_interval(2);
                   counter_window = counter_window - window_step;
                end
            end
        end
    end
    
    delete temp.mat
end
