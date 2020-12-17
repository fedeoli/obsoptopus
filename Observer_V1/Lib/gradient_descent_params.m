%% gradient descent method
function x_opt = gradient_descent_params(x_start,alpha,Niter,thresh)

    global DynOpt params
    
    % init
    iter = 1;
    
    % first step
    alpha = alpha_dynamic(x_start.b0);

    % set derivative buffer
    windowsamples = floor(DynOpt.WindowSamples/(DynOpt.Nts+1));
    checkpoint = max(0,floor(DynOpt.ActualTimeIndex/(DynOpt.Nts+1))-windowsamples);   
    if checkpoint > DynOpt.d1_derivative
        checkpoint = DynOpt.d1_derivative;
    end 
    if checkpoint > 0
       nzeros = DynOpt.d1_derivative-checkpoint;
       range = DynOpt.BackTimeIndex-DynOpt.Nts*checkpoint:DynOpt.Nts:DynOpt.BackTimeIndex-DynOpt.Nts;
       DynOpt.buf_dyhat_grad = [zeros(1,nzeros), DynOpt.OptXstory(DynOpt.params.observed_state,range)];
    end
       
    % grad computation
    [grad, xend] = sensitivity_equations_params_v3(params,x_start);
    
    % store grad
    DynOpt.grad_story(:,DynOpt.ActualTimeIndex) = grad;
    
    if norm(grad) >= thresh
        x_opt = x_start.b0 - alpha.*grad;
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
        
        
        
        % grad computation
        [grad, xend] = sensitivity_equations_params_v3(params,x_start);
        x_opt = x_opt - alpha.*grad;  
        iter = iter+1;
    end
end