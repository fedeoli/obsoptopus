%% symbolic analysis of the cost function

global DynOpt

% symbolic definition
syms wx wy wz;          % System angular velocity
syms q1 q2 q3 q0;       % Quaternion
syms n_iter

%optimization vector  
if DynOpt.bias_dyn == 0
    X = [x;DynOpt.OptXstory(end-(length(DynOpt.param_estimate)-1):end,DynOpt.BackTimeIndex)]; 
else
    X = x;
end

if DynOpt.integration_pos == 1
    X = [DynOpt.OptXstory(1:6,DynOpt.BackTimeIndex); X];
end



if n_item == length(DynOpt.Y)
    n_iter = n_item;
else
    n_iter = n_item;
end

for j=1:n_iter
    
    %evaluate the weighted in time cost function at this iteration time
    if(DynOpt.ForwardOptimization ~= 1) %backward        
        %%%%%% TO BE DONE %%%%%%
    else
        
        % get measure
        [DynOpt.buf_dyhat_temp, Yhat] = DynOpt.get_measure(X,j,1,DynOpt.buf_dyhat_temp,DynOpt.Yhat_full_story,params);
                
        J_meas = zeros(1,size(Yhat,1));
        J_der = zeros(1,size(Yhat,1));
        J_int = zeros(1,size(Yhat,1));
        J_dyn = zeros(1,size(Yhat,1));
        
        % J meas
        for i=1:length(J_meas)
            diff = (DynOpt.Y(i,DynOpt.w-n_item+j)-Yhat(i,1));
            J_meas(i) = DynOpt.scale_factor(j,1,i)*(diff)^2;
        end
        
        % J derivative
        n_der = DynOpt.dim_out;
        for i=1:n_der
            diff = (DynOpt.dY(i,DynOpt.w-n_item+j)-Yhat(i,2));
            J_der(i) = DynOpt.scale_factor(j,2,i)*(diff)^2;
        end
        
        % integral
        n_int = DynOpt.dim_out;
        for i=1:n_int
            diff = (DynOpt.intY(i,DynOpt.w-n_item+j)-Yhat(i,3));
            J_int(i) = DynOpt.scale_factor(j,3,i)*(diff)^2;
        end
        
        % J dynamics
        X_meas = X;
        Yhat_meas = get_measure_dyn_v1(X_meas,j,1,params);
        n_int = 3;
        for i=1:n_int
            diff = (DynOpt.dY(i,DynOpt.w-n_item+j)-Yhat_meas(4+i));
            J_dyn(i) = DynOpt.scale_factor(j,4,i)*(diff)^2;
        end
        
        % normalised quaternion
        quatnorm_val = quatnorm(X(1:4)');
        J_quat = DynOpt.scale_factor(j,5,i)*(1-quatnorm_val)^2;
        
        %%% temp term %%%
%         if X(1) < 0
%             J_temp = 10;
%         else
%             J_temp = 0;
%         end
        
        Jtot = Jtot + sum(J_meas) + sum(J_der) + sum(J_int) + sum(J_dyn) + J_quat;% + J_temp;     
        
        % save cunks of J
        if ~isnan(sum(J_meas)) && ~isinf(sum(J_meas))
            DynOpt.J_meas_buf = J_meas;
        end
        if ~isnan(sum(J_der)) && ~isinf(sum(J_der))
            DynOpt.J_der_buf = J_der;
        end
        if ~isnan(sum(J_int)) && ~isinf(sum(J_int))
            DynOpt.J_int_buf = J_int;
        end
        if ~isnan(sum(J_dyn)) && ~isinf(sum(J_dyn))
            DynOpt.J_dyn_buf = J_dyn;
        end
        if ~isnan(sum(J_int)) && ~isinf(sum(J_int))
            DynOpt.J_quat_buf = J_quat;
        end
                       
    end
    
end
    if n_iter > 0
        DynOpt.Yhat_temp = Yhat;
    else
        Jtot = 1; 
    end