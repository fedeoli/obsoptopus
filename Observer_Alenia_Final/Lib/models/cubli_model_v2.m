%% Inverse pendulum with flywheeò
function x_dot = cubli_model_v2(t,x,params)

    x_dot = zeros(length(x),1);
    
    x_dot(1) = x(3);
    x_dot(2) = x(4);
    
    x_dot(3) = (params.Lt*params.M*params.g*sin(x(1)) - params.u + params.Fw*x(4) - params.Fc*x(3))/params.If;
    x_dot(4) = (params.u*(params.If + params.Iw) - params.Fw*x(4)*(params.If + params.Iw) - params.Lt*params.M*params.g*sin(x(1))*params.Iw + params.Fc*x(3)*params.Iw)/(params.If*params.Iw);
end