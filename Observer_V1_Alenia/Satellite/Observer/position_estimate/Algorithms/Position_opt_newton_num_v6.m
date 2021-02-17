%% Newton Rhapson method - Symbolic analysis and numerical testing
%% Numeric analysis
% This method compute the position estimate for the entire  fleet. The
% computation is therefore centralized. 
% Input: 
%   1) Chi: (nagent x 3) matrix containing the a priori estimate for all
%   the agents.
%   2) GPS: (nagent x 3) matrix containing the GPS measures for all the
%   agents.
%   3) adjmat_UWB: (nagent x nagent) matrix ontaining the set of relative
%   distances measured by UWB. 
%   4) weights: array of parameters. It is built to be [W_UWB W_GPS W_SIGMA]
%   5) expval: array of exponents for the cost function. By default it
%   should be set to be [2 2 2].
%   6) Check_dist: boolean flag used to trigger the security check on the
%   final estimate. 
% Output:
%   1) opt: data structure containing the set of nagent position estimates
%   as well as the cost function, the gradient and hessian.
function opt = Position_opt_newton_num_v6(Chi, GPS, adjmat_UWB, weights, expval,check_dist)

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
C = zeros(1,nagent);
Chi_opt = zeros(nagent,3);
n_iter = zeros(1,nagent);
grad_C = zeros(nagent, 3);
hessian_C = zeros(3,3,nagent);

%% Init gradient and hessian - first optimization step
if get_elapsed
    tic
end

% Creation of the cost function, its gradient and its hessian - 1st iter 
for i = 1:nagent
    out = Cj_num_v7(x_hat0,Chi,GPS,adjmat_UWB,i,nagent,weights,expval);
    C(i) = out.C;
    grad_C(i,:) = out.grad;
    hessian_C(:,:,i) = out.hessian;
end
if get_elapsed
    temp = toc;
    fprintf('Function creation elapsed in %s sec\n',toc);
end

%% optimization procedure
if get_elapsed
    tic 
end

% Newton-Raphson method implementation. 
for i = 1:nagent
    matherror_flag = 0;
%     Chi_opt(i,:) = Chi(i,:);
    while (norm(grad_C(i,:)) >= thresh && n_iter(i) < max_iter)
        n_iter(i) = n_iter(i) + 1;
        
        % update step
        Chi_opt(i,:) = Chi_opt(i,:) - alpha*grad_C(i,:)*pinv(hessian_C(:,:,i));
        
        % Creation of the cost function, its gradient and its hessian. 
        out = Cj_num_v7(Chi_opt,Chi,GPS,adjmat_UWB,i,nagent,weights,expval);
        C(i) = out.C;
        grad_C(i,:) = out.grad;
        hessian_C(:,:,i) = out.hessian;
        
        % Check if Inf or NaN occurred
        if isnan(norm(grad_C(i,:))) || isinf(norm(grad_C(i,:))) 
           matherror_flag = 1;
           break
        end      
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
    for i = 1:nagent
       if  (norm(Chi_start(i,:)-Chi_opt(i,:)) >= dist_thresh)
           Chi_opt(i,:) = Chi_start(i,:);
           C(i) = -1;
           grad_C(i,:) = -1*ones(1,3);
           hessian_C(:,:,i) = -1*ones(3,3);
       end
    end
end

%% Final computations and structure assignment
grad_norm = zeros(1,nagent);
% Norm of the gradient
for i = 1:nagent
    grad_norm(i) = norm(grad_C(i,:));
end

opt.C = C;
opt.grad = grad_C;
opt.grad_norm = grad_norm;
opt.hessian = hessian_C;
opt.Chi_opt = Chi_opt;
opt.n_iter = n_iter;

end