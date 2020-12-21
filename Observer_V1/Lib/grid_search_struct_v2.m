%% Montecarlo simulations 
function grid_search_struct_v2(params_interval,pathname)
    % save data in path
    path = pathname;
    
    A_gamma = params_interval(1,1):params_interval(1,3):params_interval(1,2);
    A_gamma1 = params_interval(2,1):params_interval(2,3):params_interval(2,2);
    A_Wt = params_interval(3,1):params_interval(3,3):params_interval(3,2);
    
    B = allcomb(A_gamma,A_gamma1,A_Wt);
    n_comb = length(B);

    % grid search
    for i=1:n_comb
        file = strcat('/simulation_',int2str(i));

        final_path = strcat(path,file);
        filename = strcat(final_path,'.mat');
        save temp
        
        % launch algorithm
        keep filename path final_path i B
        init_struct
        
        struct.gamma = B(i,1);
        struct.gamma1 = B(i,2);
        struct.Wt = B(i,3);
        
        [DynOpt, params] = MainOpt_DEZ_general_v15_fun_params(struct);
        save(final_path);
        
        load temp
    end
    
    delete temp.mat
end
