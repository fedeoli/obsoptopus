%% gradient descent method
function [x_opt, J] = gradient_descent_v2(x_start,params)

    global DynOpt 
    
    % init
    iter = 1;
    thresh = DynOpt.grad_thresh;
    Niter = DynOpt.max_iter;
    
    % first step
    grad = grad_analytic(params,x_start);
    
    % select alpha
    alpha = 1e-2*norm(x_start);
    
    % first step - init solution
    if norm(grad) >= thresh
        x_opt = x_start - alpha.*grad;
    else
        x_opt = x_start;
    end
    
    % while loop
    while iter < Niter && norm(grad) >= thresh
        % recast initial condition
        x_start = x_opt;
        
        % gradient descent
        grad = grad_analytic(params,x_start);
        alpha = 1e-2*norm(x_start);
        x_opt = x_opt - alpha.*grad;     
        iter = iter+1;
    end
    
    % cost function
    J = DynOpt.cost_function(x_opt);
end