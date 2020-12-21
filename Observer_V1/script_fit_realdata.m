
%% simulate
if 1
    clear all
    pathname = 'simulation/runaway/data';
    
    params_interval = [0.2 1 0.1;...
                        1 3 0.5;...
                        0.5 5 0.5];
    grid_search_struct_v2(params_interval,pathname)
end


