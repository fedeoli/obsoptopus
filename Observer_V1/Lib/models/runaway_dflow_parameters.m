%% derivative of the flow - Khalil chapter 3
function csi2_dot = runaway_dflow_parameters(t,csi1,csi2,params)

    % dynamic matrix
    dflow = [-2*csi1(1), -2*csi1(2), 0, 0;...
             params.gamma*csi1(2), -params.ni + params.gamma*csi1(1)-params.gamma1/(1+csi1(2)/params.Wt)+params.gamma1*csi1(2)/(params.Wt*(1+csi1(2)/params.Wt)^2),...
                                   csi1(1)*csi1(2)+params.S, -csi1(2)/(1+csi1(2)/params.Wt);...
             0, 0, 0, 0;...
             0, 0, 0, 0];
         
    % consider the eps_coef
    dflow = params.eps_coef*dflow;
    
    % evolution
    csi2_dot = dflow*csi2;
end