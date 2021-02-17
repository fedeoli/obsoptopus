%% gradient descent method
function x_opt = gradient_descent_params_v4(x_start,Niter,thresh)

    global DynOpt params
    
    % init
    iter = 1;
    
    % read measure 
    measure_forward = 1;
    
    % first step
    alpha = alpha_dynamic(x_start.b0);

    % set derivative buffer
    gradstory = DynOpt.OptXstory;
       
    % grad computation
    [grad, xend] = sensitivity_equations_params_v5(params,x_start);
    
    % store grad
    DynOpt.grad_story(:,end+1) = grad;
    
    if norm(grad) >= thresh
        x_opt = x_start.b0 - alpha.*grad;
        
        % store and propagate
        % params and state update
        if DynOpt.identify == 1
            params_update(x_opt);
        end
        
        % update 
        x_propagate = x_opt;

        state_update;
    else
        x_opt = x_start.b0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% 
    % GRADIENT DESCENT WHILE LOOP 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    while iter < Niter && norm(grad) >= thresh
        % recast initial condition
        xstart.b0 = x_opt;
        xstart.e0 = xend.e;
        
        % gradient descent
        alpha = alpha_dynamic(x_start.b0);
        params_update(x_opt);
        
        % grad computation
        [grad, xend] = sensitivity_equations_params_v5(params,x_start);
        x_opt = x_opt - alpha.*grad;  
        
        % store grad
        DynOpt.grad_story(:,end+1) = grad;
        
        % store and propagate
        % params and state update
        if DynOpt.identify == 1
            params_update(x_opt);
        end
        
        % update 
        state_update;
                
        % iter
        iter = iter+1;
    end
end