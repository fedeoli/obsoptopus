%% Chapter 3 Khalil - sensitivity equations
function [J_grad, xend] = sensitivity_equations_params_v5(params,x0)

    global DynOpt 
    
    J_grad = zeros(1,DynOpt.StateDim+DynOpt.nparams);
    
    % set the derivative buffer as before the optimisation process (multiple f computation)
    if DynOpt.BackTimeIndex >= DynOpt.d1_derivative
        DynOpt.buf_dyhat_grad = DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.d1_derivative-1):DynOpt.BackTimeIndex);
    else
        init_pos = DynOpt.d1_derivative-DynOpt.BackTimeIndex;
        DynOpt.buf_dyhat_grad = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.BackTimeIndex-1):DynOpt.BackTimeIndex)];
    end
    if DynOpt.BackTimeIndex >= DynOpt.d1_dderivative
        DynOpt.buf_ddyhat_grad = DynOpt.dYhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.d1_dderivative-1):DynOpt.BackTimeIndex);
    else
        init_pos = DynOpt.d1_dderivative-DynOpt.BackTimeIndex;
        DynOpt.buf_ddyhat_grad = [zeros(DynOpt.dim_out,init_pos), DynOpt.dYhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.BackTimeIndex-1):DynOpt.BackTimeIndex)];
    end
    
    for i=1:(DynOpt.w)       
       % measured output
       a(1,1) = DynOpt.Y(i);
       a(2,1) = DynOpt.dY(i);
       a(3,1) = DynOpt.ddY(i);
       
       % propagate the flow and the Jacobian of the flow
       b = x0.b0;
       e = x0.e0;
       
       for j=1:DynOpt.Nts*(i-1)
          % state propagation
          b = PlantJumpMap_general_notime_params(b,DynOpt.model,DynOpt.ForwardOptimization,params);
          e = PlantJumpMap_general_notime_params_hybrid(b,e,DynOpt.model_flow,DynOpt.ForwardOptimization,params);
          
          % derivative propagation
          % map the flow on the output
          c = b(params.observed_state);
          [DynOpt.buf_dyhat_grad, dc] = IterativePseudoDerivative(DynOpt.Ts,c,DynOpt.c1_derivative,DynOpt.d1_derivative,0,DynOpt.buf_dyhat_grad);
          [DynOpt.buf_ddyhat_grad, ddc] = IterativePseudoDerivative(DynOpt.Ts,dc,DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,DynOpt.buf_ddyhat_grad);
          c_tot = [c; dc; ddc];
       end
       
       if i==1
           % derivative propagation
           % map the flow on the output
           c = b(params.observed_state);
           [DynOpt.buf_dyhat_grad, dc] = IterativePseudoDerivative(DynOpt.Ts,c,DynOpt.c1_derivative,DynOpt.d1_derivative,0,DynOpt.buf_dyhat_grad);
           [DynOpt.buf_ddyhat_grad, ddc] = IterativePseudoDerivative(DynOpt.Ts,dc,DynOpt.c1_dderivative,DynOpt.d1_dderivative,0,DynOpt.buf_ddyhat_grad);
           c_tot = [c; dc; ddc];
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % gradient of the output mapping h(x)
       % this is analytic so change it manually
       d = zeros(DynOpt.y_end,DynOpt.nparams+DynOpt.StateDim);
       % first row [0 1 0 0]
       d(1,params.observed_state) = 1;
       % second row
       d(2,1) = params.gamma*b(2);
       d(2,2) = -params.ni+params.gamma*b(1)-params.gamma1*(params.Wt^2+params.Wt*b(2)-params.Wt^2*b(2))/(params.Wt+b(2))^2;
       d(2,3) = b(1)*b(2)+params.S;
       d(2,4) = -(params.Wt*b(2))/(params.Wt+b(2));
       % third row
       d(3,1) = -params.ni*d(2,1)-2*params.gamma*b(2)^2+params.gamma*b(1)*d(2,1)-params.gamma1*d(2,1)*(params.Wt^2+params.Wt*b(2)-params.Wt^2*b(2))/(params.Wt+b(2))^2;
       d(3,2) = -params.ni*d(2,2)-params.gamma*b(2)*(2*b(1))-params.gamma*b(3)+params.gamma1*b(1)*d(2,2)-params.gamma1*d(2,2)*(params.Wt^2+params.Wt*b(2)-params.Wt^2*b(2))/(params.Wt+b(2))^2 ...
                -params.gamma1*b(4)*((params.Wt-params.Wt^2)*(params.Wt+b(2))^2 - (params.Wt^2+params.Wt*b(2)-params.Wt^2*b(2))*2*(params.Wt+b(2)))/(params.Wt+b(2))^4;
       d(3,3) = b(3)*b(2) + b(1)*b(4);
       d(3,4) = -b(4)*(params.Wt^2+params.Wt*b(2)-params.Wt^2*b(2))/(params.Wt+b(2))^2;
       
       % add the eps_coef
       d = params.eps_coef*d;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
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