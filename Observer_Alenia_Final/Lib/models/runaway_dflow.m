%% derivative of the flow - Khalil chapter 3
function csi2_dot = runaway_dflow(t,csi1,csi2,params)

    % dynamic matrix
    dflow = [-2*csi1(1), -2*csi1(2);...
             params.gamma*csi1(2), -params.ni + params.gamma*csi1(1)-params.gamma1/(1+csi1(2)/params.Wt)+params.gamma1*csi1(2)/(params.Wt*(1+csi1(2)/params.Wt)^2)];
    % evolution
    csi2_dot = dflow*csi2;
end