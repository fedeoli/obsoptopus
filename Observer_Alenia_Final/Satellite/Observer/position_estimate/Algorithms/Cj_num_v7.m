%% Cost function creation - Newton Raphson method
% This method computes cost function, gradient and hessian assuming to
% completely know Chi,GPS,UWB of all agents. This method can be used to
% implement a centralized simulation of the fleet.
% Input: 
%   1) x_hat: (nagent x 3) matrix with the current position estimates for
%   all the agents.
%   2) Chi: (nagent x 3) matrix containing the a priori estimate for all
%   the agents.
%   3) GPS: (nagent x 3) matrix containing the GPS measures for all the
%   agents.
%   4) adjmat_UWB: (nagent x nagent) matrix ontaining the set of relative
%   distances measured by UWB. 
%   5) j_select: number of the agent considered.
%   6) weights: array of parameters. It is built to be [W_UWB W_GPS W_SIGMA]
%   7) expval: array of exponents for the cost function. By default it
%   should be set to be [2 2 2].
%   8) Check_dist: boolean flag used to trigger the security check on the
%   final estimate. 
% Output:
%   1) opt: data structure containing the set of nagent position estimates
%   as well as the cost function, the gradient and hessian.
function out = Cj_num_v7(x_hat, Chi, GPS, adjmat, j_select, nagent, weights, expval)
%%  extract data from input 
%   extract weights from input
    Wuwb = weights(1); 
    Wgps = weights(2);
    Wsigma = weights(3);
    
%   extract exponents from input
    rho_uwb = expval(1);
    rho_gps = expval(2);
    rho_sigma = expval(3);
    
%% CREATE COST FUNCTION
%  create first term J1
    Juwb = 0;
    for i=1:nagent
        Juwb = Juwb + ( norm(x_hat(j_select,:) - Chi(i,:)) - adjmat(i,j_select) )^rho_uwb;
    end
    Juwb = Wuwb*Juwb;

%  create second term J2
   Jgps = Wgps*norm(x_hat(j_select,:)-GPS(j_select,:))^rho_gps;

%  create third term J3
   Jsigma = Wsigma*norm(x_hat(j_select,:)-Chi(j_select,:))^rho_sigma;

%  Assign Cj (cost function)
   Cj = Juwb + Jgps + Jsigma;
    
%% COMPUTE GRADIENT
grad_Juwb = zeros(1,3);
grad_Jgps = zeros(1,3);
grad_Jsigma = zeros(1,3);
for j = 1:3
   for i = 1:nagent
       T = (norm(x_hat(j_select,:) - Chi(i,:)) - adjmat(i,j_select))*(x_hat(j_select,j)-Chi(i,j))/norm(x_hat(j_select,:) - Chi(i,:));
       if ~isnan(T)
            grad_Juwb(j) = grad_Juwb(j) + T;
       end
   end
   grad_Juwb(j) = 2*Wuwb*grad_Juwb(j);
   grad_Jgps(j) = 2*Wgps*(x_hat(j_select,j)-GPS(j_select,j));
   grad_Jsigma(j) = 2*Wsigma*(x_hat(j_select,j)-Chi(j_select,j));
end
grad_Cj = grad_Juwb + grad_Jgps + grad_Jsigma;

%% COMPUTE HESSIAN
hessian_Juwb = zeros(3,3);
hessian_Jgps = zeros(3,3);
hessian_Jsigma = zeros(3,3);
for h = 1:3
   % HESSIAN JUWB
   for k = 1:3
      for i = 1:nagent
        % first term
        T11 = (x_hat(j_select,h)-Chi(i,h))*(x_hat(j_select,k)-Chi(i,k))/norm(x_hat(j_select,:)-Chi(i,:));
        if h == k
            T12 =  norm(x_hat(j_select,:)-Chi(i,:)) - adjmat(i,j_select);
        else
            T12 = 0;
        end
        T1 = norm(x_hat(j_select,:)-Chi(i,:))*(T11 + T12);     
        % second term
        T2 = (norm(x_hat(j_select,:)-Chi(i,:)) - adjmat(i,j_select))*(x_hat(j_select,h)-Chi(i,h))*(x_hat(j_select,k)-Chi(i,k))/norm(x_hat(j_select,:)-Chi(i,:));        
        % combine two terms
        T = (T1 - T2)/(norm(x_hat(j_select,:)-Chi(i,:)))^2;
        if ~isnan(T)
            hessian_Juwb(h,k) = hessian_Juwb(h,k) + T;
        end
      end
      hessian_Juwb(h,k) = 2*Wuwb*hessian_Juwb(h,k);
      
      % HESSIAN JGPS + JSIGMA
      if h == k
        hessian_Jgps(h,k) = 2*Wgps;
        hessian_Jsigma(h,k) = 2*Wsigma;
      end
      
   end
end
hessian_Cj = hessian_Juwb + hessian_Jgps + hessian_Jsigma;

out.C = Cj;
out.grad = grad_Cj;
out.hessian = hessian_Cj;
    
end