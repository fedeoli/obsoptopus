%% MONTECARLO SIMULATION
function montecarlo_run(n_sim)

    % cycle over the simulations
    for i=1:n_sim
        
       % init the simulation
       init_struct;
       
       % randomly define the parameters
       struct.params_init.gamma = (struct.ub_init(3)-struct.lb_init(3)).*rand(1) + struct.lb_init(3);
       struct.params_init.gamma1 = (struct.ub_init(4)-struct.lb_init(4)).*rand(1) + struct.lb_init(4);
       
       % call optimiser
       [DynOpt, params] = MainOpt_DEZ_general_v22_fun_params(struct);
       
       file = strcat('simulation/montecarlo/simulation_',int2str(i));
       save(file);
       
    end

end