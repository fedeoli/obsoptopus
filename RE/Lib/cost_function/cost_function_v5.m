%% cost function
function Jtot = cost_function_v5(x,params)

global DynOpt

%x: the estimated optimal state and parameters vector provided by fmin...

%DynOpt.Ts : Sampling Time
%DynOpt.Nts : Inter-sampling Time
%DynOpt.X : actual state (t_k), a column vector, first the state, then stacked the parameters
%DynOpt.U : input vector, last value at time t_k, a matrix with w column vectors with the number of rows as the number of "instantaneous" input
%DynOpt.Y : output vector, last value at time t_k, a matrix with w column vectors with the number of rows as the number of "instantaneous" outputs
%DynOpt.w : window size
%DynOpt.StateDim : number of state variables (subject to dynamics, not parameters)
%DynOpt.Weight : window weight size, last on the more recent t_k

% cost function init
Jtot = 0;

% optimise on the euler angles
temp = (RotationConversion_V2_1('EA321toQ', DynOpt.temp_x0(1:3)'))';
x = [temp; x(4:end)];

%optimization vector  
if DynOpt.bias_dyn == 0
    X = [x;DynOpt.OptXstory(end-(length(DynOpt.param_estimate)-1):end,DynOpt.BackTimeIndex)]; 
else
    X = x;
end

% params update
if DynOpt.identify == 1
    params = DynOpt.params_update(DynOpt.X,params);
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
            J_meas(i) = DynOpt.scale_factor(1,i)*(diff)^2;
        end
        
        % J derivative
        n_der = DynOpt.dim_out;
        for i=1:n_der
            diff = (DynOpt.dY(i,DynOpt.w-n_item+j)-Yhat(i,2));
            J_der(i) = DynOpt.scale_factor(2,i)*(diff)^2;
        end
        
        % integral
        n_int = DynOpt.dim_out;
        for i=1:n_int
            diff = (DynOpt.intY(i,DynOpt.w-n_item+j)-Yhat(i,3));
            J_int(i) = DynOpt.scale_factor(3,i)*(diff)^2;
        end
        
        % J dynamics
        X_meas = X;
        Yhat_meas = get_measure_dyn_v1(X_meas,j,1,params);
        n_int = DynOpt.dim_out;
        for i=1:n_int
            diff = (DynOpt.dY(i,DynOpt.w-n_item+j)-Yhat_meas(4+i));
            J_dyn(i) = DynOpt.scale_factor(4,i)*(diff)^2;
        end
        
        % normalised quaternion
        quatnorm_val = quatnorm(X(1:4)');
        J_quat = DynOpt.scale_factor(5,1)*(1-quatnorm_val)^2;
        
        Jtot = Jtot + sum(J_meas) + sum(J_der) + sum(J_int) + sum(J_dyn) + J_quat;     
        
        % save J_meas and J_der
        if ~isnan(sum(J_meas)) && ~isinf(sum(J_meas))
            DynOpt.J_meas(end+1) = sum(J_meas);
        end
        if ~isnan(sum(J_der)) && ~isinf(sum(J_der))
            DynOpt.J_der(end+1) = sum(J_der);
        end
        if ~isnan(sum(J_int)) && ~isinf(sum(J_int))
            DynOpt.J_int(end+1) = sum(J_int);
        end
                       
    end
    
end
    if n_iter > 0
        DynOpt.Yhat_temp = Yhat;
    else
        Jtot = 1; 
    end
end
    