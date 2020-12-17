%% Chapter 3 Khalil - sensitivity equations
function [J_grad, xend] = sensitivity_equations_params_v3(params,x0)

    global DynOpt 
    
    J_grad = zeros(1,DynOpt.StateDim+DynOpt.nparams);
    
    % set the derivative buffer as before the optimisation process (multiple f computation)
    if DynOpt.BackTimeIndex >= DynOpt.d1_derivative
        DynOpt.buf_dyhat_grad = DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.d1_derivative-1):DynOpt.BackTimeIndex);
    else
        init_pos = DynOpt.d1_derivative-DynOpt.BackTimeIndex;
        DynOpt.buf_dyhat_grad = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.BackTimeIndex-1):DynOpt.BackTimeIndex)];
    end
    
    for i=1:(DynOpt.w)       
       % measured output
       a = DynOpt.Y(1:DynOpt.y_end,i);
       
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
       d = zeros(1,DynOpt.nparams+DynOpt.StateDim);
       d(1,params.observed_state) = 1;
                
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