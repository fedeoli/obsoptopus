
%% simulate
if 0
    clear all
    pathname = 'simulation/runaway/data_noopt';
    
    params_interval = [0.2 1.5 0.1;...
                        1 5 0.5;...
                        0.5 5 0.5];
    grid_search_struct_v2(params_interval,pathname)
end

%% analyse
if 0
    pathname = 'simulation/runaway/data_noopt';
    
    params_interval = [0.2 1.5 0.1;...
                        1 5 0.5;...
                        0.5 5 0.5];
    
    plot_flag = 1;

    out = grid_analyse_v2(params_interval,pathname);
    
end

%% multiple plot
if 1
    start = 1253;
    stop = 1255;

    figure
    hold on
    for i=start:stop
        plot(out.state(i,:));
    end
end

