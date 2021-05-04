%% gradient descent method
function x_opt = gradient_descent_params_v2(x_start,alpha,Niter,thresh)

    global DynOpt params
    
    % init
    iter = 1;
    
    % first step
    alpha = alpha_dynamic(x_start.b0);

    % set derivative buffer
    gradstory = DynOpt.OptXstory;
    
    % windowsamples
    windowsamples = floor(DynOpt.WindowSamples/(DynOpt.Nts+1));   
    checkpoint = max(0,floor(DynOpt.ActualTimeIndex/(DynOpt.Nts+1))-windowsamples);
    
    if checkpoint > DynOpt.d1_derivative
        checkpoint = DynOpt.d1_derivative;
    end 
    if checkpoint > 0
       nzeros = DynOpt.d1_derivative-checkpoint;
       range = DynOpt.BackTimeIndex-DynOpt.Nts*checkpoint:DynOpt.Nts:DynOpt.BackTimeIndex-DynOpt.Nts;
       DynOpt.buf_dyhat_grad = [zeros(1,nzeros), gradstory(DynOpt.params.observed_state,range)];
    end
       
    % grad computation
    [grad, xend] = sensitivity_equations_params_v3(params,x_start);
    
    % store grad
    DynOpt.grad_story(:,end+1) = grad;
    
    if norm(grad) >= thresh
        x_opt = x_start.b0 - alpha.*grad;
        
        % update gradstory
        params_update(x_opt);
        x_propagate = x_opt;
        for j =1:DynOpt.WindowSamples-1     %in absolute time as: k-(w-1)*Nts+1:k,
            set_input(DynOpt.BackTimeIndex+j);
            x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,1,params);
            gradstory(:,DynOpt.BackTimeIndex+j) = x_propagate;
        end
    else
        x_opt = x_start.b0;
    end
    
    % while loop
    while iter < Niter && norm(grad) >= thresh
        % recast initial condition
        xstart.b0 = x_opt;
        xstart.e0 = xend.e;
        
        % gradient descent
        alpha = alpha_dynamic(x_start.b0);
        params_update(x_opt);
        
        % set derivative buffer
        windowsamples = floor(DynOpt.WindowSamples/(DynOpt.Nts+1));
        checkpoint = max(0,floor(DynOpt.ActualTimeIndex/(DynOpt.Nts+1))-windowsamples);  
        
        if checkpoint > DynOpt.d1_derivative
            checkpoint = DynOpt.d1_derivative;
        end 
        if checkpoint > 0
           nzeros = DynOpt.d1_derivative-checkpoint;
           range = DynOpt.BackTimeIndex-DynOpt.Nts*checkpoint:DynOpt.Nts:DynOpt.BackTimeIndex-DynOpt.Nts;
           DynOpt.buf_dyhat_grad = [zeros(1,nzeros), gradstory(DynOpt.params.observed_state,range)];
        end
        
        % grad computation
        [grad, xend] = sensitivity_equations_params_v3(params,x_start);
        x_opt = x_opt - alpha.*grad;  
        
        % update gradstory
        params_update(x_opt);
        x_propagate = x_opt;
        for j =1:DynOpt.WindowSamples-1     %in absolute time as: k-(w-1)*Nts+1:k,
            set_input(DynOpt.BackTimeIndex+j);
            x_propagate =  PlantJumpMap_general_notime_params(x_propagate,DynOpt.model,1,params);
            gradstory(:,DynOpt.BackTimeIndex+j) = x_propagate;
        end
        
        % iter
        iter = iter+1;
    end
end