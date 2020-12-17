%% gradient descent method
function x_opt = gradient_descent(x_start,alpha,Niter,thresh,params)

    global DynOpt
    
    % init
    iter = 1;
    
    % first step
    [grad, xend] = sensitivity_equations(params,x_start);
    
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
        [grad, xend] = sensitivity_equations(params,x_start);
        x_opt = x_opt - alpha.*grad;     
        iter = iter+1;
    end
end