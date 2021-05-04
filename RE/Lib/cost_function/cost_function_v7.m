%% cost function
function Jtot = cost_function_v7(x,params)

% global vars
global DynOpt

% cost function init
Jtot = 0;

%optimization vector  
if DynOpt.nparams == 0
    X = [x;DynOpt.OptXstory(end-(length(DynOpt.param_estimate)-1):end,DynOpt.BackTimeIndex)]; 
else
    X = x;
end

if DynOpt.integration_pos == 1
    X = [DynOpt.OptXstory(1:6,DynOpt.BackTimeIndex); X];
end
offset = DynOpt.integration_pos*6;

% params update
if DynOpt.identify == 1
    params = DynOpt.params_update(X,params);
end

% set the derivative buffer as before the optimisation process (multiple f computation)
if DynOpt.BackTimeIndex >= DynOpt.d1_derivative
    DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.d1_derivative-1):DynOpt.BackTimeIndex);
else
    init_pos = DynOpt.d1_derivative-DynOpt.BackTimeIndex;
    DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,DynOpt.BackTimeIndex-(DynOpt.BackTimeIndex-1):DynOpt.BackTimeIndex)];
end

n_item = length(find(min(abs(DynOpt.Y),[],1)));

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
        X_meas = X();
        Yhat_meas = get_measure_dyn_v1(X_meas,j,1,params);
        n_int = 3;
        for i=1:n_int
            diff = (DynOpt.dY(i,DynOpt.w-n_item+j)-Yhat_meas(4+i));
            J_dyn(i) = DynOpt.scale_factor(j,4,i)*(diff)^2;
        end
        
        % normalised quaternion
        quatnorm_val = quatnorm(X(offset+1:offset+4)');
        J_quat = DynOpt.scale_factor(j,5,i)*(1-quatnorm_val)^2;
        
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
end
    