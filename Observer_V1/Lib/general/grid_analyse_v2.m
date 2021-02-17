%% Montecarlo simulations 
function struct = grid_analyse_v2(params_interval,pathname)
    
    %% Data analysis
    path = pathname;
    
    A_gamma = params_interval(1,1):params_interval(1,3):params_interval(1,2);
    A_gamma1 = params_interval(2,1):params_interval(2,3):params_interval(2,2);
    A_Wt = params_interval(3,1):params_interval(3,3):params_interval(3,2);
    
    B = allcomb(A_gamma,A_gamma1,A_Wt);
    n_comb = length(B);
    
    mean_error = -1;
        
    % grid search
    for i=1:n_comb

        file = strcat('/simulation_',int2str(i));

        final_path = strcat(path,file);
        load(final_path);

        if DynOpt.ObserverOn == 1
            mean_error(end+1) = norm(sum(abs(DynOpt.OptErrorStory),2));
            keep path mean_error B 
        else
            state(i,:) = DynOpt.stateStory(2,:);
            keep path B state
        end
           
    end
    
    if exist('mean_error','var')
        struct.mean_error = mean_error(2:end); 
    end
    
    if exist('state','var')
        struct.state = state;
    end
    
    struct.B = B;
end