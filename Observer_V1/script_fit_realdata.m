
%% simulate
if 1
    clear all
    pathname = 'simulation/runaway/data_noopt';
    
    params_interval = [0.2 1 0.1;...
                        1 3 0.5;...
                        0.5 5 0.5];
    grid_search_struct_v2(params_interval,pathname)
end

%% analyse
if 1
    pathname = 'simulation/runaway/data_noopt';
    
    params_interval = [0.2 1 0.1;...
                        1 3 0.5;...
                        0.5 5 0.5];
    
    plot_flag = 1;

    out = grid_analyse_v2(params_interval,pathname);
    
end

%% multiple plot
if 0
    start = 220;
    stop = 250;

    figure
    hold on
    for i=1:start:stop
        plot(out.state(i,:));
    end
end


