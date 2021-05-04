%% setup contraints for opt problem
function [c,ceq] = setup_constraints(x)
    global DynOpt
    n_params = length(DynOpt.param_estimate);
    
    % inequalities
    c = [];
%     c = zeros(n_params,1);
%     for i=1:n_params    
%         c(i) = -x(end-(i-1));
%     end
    
    % equalities
    ceq = [];
end