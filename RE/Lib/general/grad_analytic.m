%% Chapter 3 Khalil - sensitivity equations
function J_grad = grad_analytic(params,x0)
    global DynOpt
    
    state = 1e-1*x0;
    J_grad = state;
end