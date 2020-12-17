%% Chapter 3 Khalil - sensitivity equations
function [J_grad, xend] = sensitivity_equations(params,x0)

    global DynOpt 
    
    % grad init
    n_obs = length(DynOpt.params.observed_state);
    J_grad = zeros(1,n_obs);
    
    % grad computation
    for i=1:DynOpt.w       
       % measured output
       a = DynOpt.Y(1,i);
       
       % propagate the flow and the Jacobian of the flow
       b = x0.b0;
       e = x0.e0;
       for j=1:DynOpt.Nts*i
          b = PlantJumpMap_general_notime_params(b,DynOpt.model,DynOpt.ForwardOptimization,params);
          e = PlantJumpMap_general_notime_params_hybrid(b,e,DynOpt.model_flow,DynOpt.ForwardOptimization,params);
       end
       
       % map the flow on the output
       c = b(params.observed_state);
       
       % gradient of the output mapping h(x)
       % this is analytic so change it manually
       d = zeros(1,DynOpt.StateDim);
       d(params.observed_state) = 1;
       
       % update J_grad
       temp = (a-c)'*d*e;
       J_grad = J_grad + temp;
    end
    
    % grad finalisation
    len_tmp = (DynOpt.StateDim+length(DynOpt.param_estimate))-length(J_grad);
    J_grad = [-J_grad, zeros(1,len_tmp)]';
    
    % store xend
    xend.b = b;
    xend.e = e;
end