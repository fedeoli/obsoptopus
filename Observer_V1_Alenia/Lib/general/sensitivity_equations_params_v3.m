%% Chapter 3 Khalil - sensitivity equations
function [J_grad, xend] = sensitivity_equations_params_v3(params,x0)

    global DynOpt 
    
     J_grad = zeros(1,DynOpt.StateDim+DynOpt.n_param_est);
    for i=1:(DynOpt.w)       
       % measured output
       a = DynOpt.Y(1:2,i);
       
       % propagate the flow and the Jacobian of the flow
       b = x0.b0;
       e = x0.e0;
       
       for j=1:DynOpt.Nts*(i-1)
          b = PlantJumpMap_general_notime_params(b,DynOpt.model,DynOpt.ForwardOptimization,params);
          e = PlantJumpMap_general_notime_params_hybrid(b,e,DynOpt.model_flow,DynOpt.ForwardOptimization,params);
       end
       
       % map the flow on the output
       c = b(params.observed_state);
       [DynOpt.buf_dyhat_grad, dc] = IterativePseudoDerivative(DynOpt.Ts,c,DynOpt.c1_derivative,DynOpt.d1_derivative,0,DynOpt.buf_dyhat_grad);
       c_tot = [c; dc];
       
       % gradient of the output mapping h(x)
       % this is analytic so change it manually
       d = zeros(2,DynOpt.StateDim+DynOpt.n_param_est);
       d(1,params.observed_state) = 1;
       d(2,:) = [params.gamma*b(2), -params.ni + params.gamma*b(1)-params.gamma1/(1+b(2)/params.Wt)+params.gamma1*b(2)/(params.Wt*(1+b(2)/params.Wt)^2),...
                                   b(1)*b(2)+params.S, -b(2)/(1+b(2)/params.Wt)];
                               
       % update J_grad
       temp = (a(1:DynOpt.y_end)-c_tot(1:DynOpt.y_end))'*d(1:DynOpt.y_end,:)*e;
       
%        if norm(temp) ~= 0
%             disp('ARARMAX');
%        end

       J_grad = J_grad + temp;
    end
    J_grad = - J_grad';
    
    % store xend
    xend.b = b;
    xend.e = e;
end