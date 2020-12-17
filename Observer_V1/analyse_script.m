
%% simulate
if 0
    pathname = 'simulation/runaway/data';
    window_interval = [5 15];
    window_step = 2;
    sample_interval = [4 16];
    sample_step = 2;

    grid_search_struct(window_interval,window_step,sample_interval,sample_step,pathname)
end

%% analyse
if 1
%     pathname = 'simulation/cubli_3s_forward';
    pathname = 'simulation/runaway/data';
    
    %%%%%% TOKAMAK %%%%%%
    window_interval = [5 15];
    window_step = 2;
    sample_interval = [4 16];
    sample_step = 2;
    
    plot_flag = 1;

    out = grid_analyse(window_interval,window_step,sample_interval,sample_step,pathname,plot_flag);
    
end



    %% select best
if 1
    % init
%     pathname = 'simulation/cubli_3s_forward';
    pathname = 'simulation/tokamak_1s_backward';
    
    %%%%%% TOKAMAK %%%%%%
    window_interval = [5 15];
    window_step = 2;
    sample_interval = [4 16];
    sample_step = 2;

    window_array = out.window_interval(1):window_step:out.window_interval(2);
    sample_array = out.sample_interval(1):sample_step:out.sample_interval(2);
    
%     window_best = window_array(out.row_lambda);
%     sample_best = sample_array(out.col_lambda);
    
    window_best = window_array(out.row_error);
    sample_best = sample_array(out.col_error);
end

%% test
if 1
    
    struct = init_struct_analyse_v2(window_best,sample_best);

    [DynOpt, params] = MainOpt_DEZ_general_v12_fun_params(struct);

end

