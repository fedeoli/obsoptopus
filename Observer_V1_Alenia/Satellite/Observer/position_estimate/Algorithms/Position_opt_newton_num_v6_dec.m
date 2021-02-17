%% Newton Rhapson method - Symbolic analysis and numerical testing
%% Numeric analysis
% This method compute the position estimate for the entire  fleet. The
% computation is therefore centralized. 
% Input: 
%   1) Chi: (nagent x 3) matrix containing the a priori estimate for all
%   the agents.
%   2) GPS: (1 x 3) matrix containing the GPS measures for the agent
%   considered.
%   3) adjmat_UWB: (nagent x nagent) matrix ontaining the set of relative
%   distances measured by UWB. 
%   4) j_select: number of the agent considered.
%   5) weights: array of parameters. It is built to be [W_UWB W_GPS W_SIGMA]
%   6) expval: array of exponents for the cost function. By default it
%   should be set to be [2 2 2].
%   7) Check_dist: boolean flag used to trigger the security check on the
%   final estimate. 
% Output:
%   1) opt: data structure containing the set of nagent position estimates
%   as well as the cost function, the gradient and hessian.
function opt = Position_opt_newton_num_v6_dec(Chi, GPS, adjmat_UWB, j_select, weights, expval,check_dist)

%% Init section
% n of agents
nagent = size(Chi,1);

% first position
Chi_start = Chi;
x_hat0 = Chi;
dist_thresh = 5e-2;

% optimization data
% Newton update - from x_k to x_k+1
thresh = 1e-4;
max_iter = 40;
alpha = 1e0;

% if this is true the algorithm shows the execution times.
get_elapsed = 0;

% accumulators
n_iter = 0;

%% Init gradient and hessian - first optimization step
if get_elapsed
    tic
end

% Creation of the cost function, its gradient and its hessian - 1st iter 
out = Cj_num_v7_dec(x_hat0,Chi,GPS,adjmat_UWB,j_select,nagent,weights,expval);
C = out.C;
grad_C = out.grad;
hessian_C = out.hessian;
if get_elapsed
    temp = toc;
    fprintf('Function creation elapsed in %s sec\n',toc);
end

%% optimization procedure
if get_elapsed
    tic 
end

% Newton-Raphson method implementation. 
matherror_flag = 0;
Chi_opt = Chi(j_select,:);
while (norm(grad_C) >= thresh && n_iter < max_iter)
    n_iter = n_iter + 1;

    % update step
    Chi_opt = Chi_opt - alpha*grad_C*pinv(hessian_C);

    % Creation of the cost function, its gradient and its hessian. 
    out = Cj_num_v7_dec(Chi_opt,Chi,GPS,adjmat_UWB,j_select,nagent,weights,expval);
    C = out.C;
    grad_C = out.grad;
    hessian_C = out.hessian;

    % Check if Inf or NaN occurred
    if isnan(norm(grad_C)) || isinf(norm(grad_C)) 
       matherror_flag = 1;
       break
    end      
end
if get_elapsed
    temp = toc;
    fprintf('Function optimization elapsed in %s sec\n',toc);
end

% alert if matherror
if matherror_flag == 1
   fprintf('Optimization blocked after matherror\n'); 
end

% check if the new position is too far from the starting point
if check_dist == 1
   if  (norm(Chi_start-Chi_opt) >= dist_thresh)
       Chi_opt = Chi_start;
       C = -1;
       grad_C = -1*ones(1,3);
       hessian_C = -1*ones(3,3);
   end
end

%% Final computations and structure assignment
% Norm of the gradient
grad_norm = norm(grad_C);

opt.C = C;
opt.grad = grad_C;
opt.grad_norm = grad_norm;
opt.hessian = hessian_C;
opt.Chi_opt = Chi_opt;
opt.n_iter = n_iter;

end