%% Runaway electrons model
% x = [T W]'
function x_dot = runaway_model(t,x,params)

    % state init
    x_dot = zeros(length(x),1);
    
    % just for reading ease
    T = x(1);
    W = x(2);

    x_dot(1) = params.eps_coef*(-2*T*W - 2*params.S + params.Q);
    x_dot(2) = params.eps_coef*(-params.ni*W + params.gamma*(T*W + params.S) - params.gamma1*W/(1+W/params.Wt));
end