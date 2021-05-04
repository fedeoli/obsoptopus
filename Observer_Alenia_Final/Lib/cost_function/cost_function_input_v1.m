%% Chapter 3 Khalil - sensitivity equations
function Jtot = cost_function_input_v1(theta_u,params)

    global DynOpt 
    
    % cost function init
    Jtot = 0;
    
    %optimization vector  
    b0 = [DynOpt.temp_x0(1:DynOpt.StateDim+DynOpt.nparams); theta_u; DynOpt.temp_x0(end)];
    X = b0;
    
    % minimum eigenvalue
    DynOpt.min_eig_past = Inf;
    
    % params update
    if DynOpt.identify == 1
        params = params_update_local(X,params);
    end
    
    n_item = length(find(min(abs(DynOpt.Y),[],1)));

    if n_item == length(DynOpt.Y)
        n_iter = n_item;
    else
        n_iter = n_item;
    end
    
    for j=1:n_iter
        
        % get adaptive buffer distance. The zero is added to make the vector of the
        % same length of the buffer as the diff command computes ((i+1)th-ith)
        buf_dist = [max(DynOpt.Y_space(1)-DynOpt.BackTimeIndex,0), diff(DynOpt.Y_space)];
        % get first nonzero element
        [~,pos] = find(buf_dist);
        % check if exist a zero
        zero_flag = pos(1)>1;
        
        %evaluate the weighted in time cost function at this iteration time
        if(DynOpt.ForwardOptimization ~= 1) %backward        
            %%%%%% TO BE DONE %%%%%%
        else
           % propagate the flow and the Jacobian of the flow
           b = b0;
           e = DynOpt.e_flow(:,:,DynOpt.BackTimeIndex);

           for i=1:sum(buf_dist(1:pos(j)))-zero_flag
               % integration
                [b, ~] = model_propagate_local(DynOpt.BackTimeIndex,DynOpt.Ts,b,params);
                e = PlantJumpMap_general_notime_params_hybrid(b,e,DynOpt.model_flow,DynOpt.ForwardOptimization,params);
                DynOpt.e_flow(:,:,DynOpt.BackTimeIndex+i) = e;
           end 
        end
        
        %%%%%%% COMPUTE JTOT HERE %%%%%%%%
        %%% update cost function %%%
        csi2 = e(5:DynOpt.StateDim,DynOpt.StateDim+1:DynOpt.StateDim+DynOpt.nparams);
        min_eig = min(eig(csi2));
%         csi2 = e(1:DynOpt.StateDim,1:DynOpt.StateDim+DynOpt.nparams);
%         min_eig = min(eig(csi2*csi2'));
        
        % state change barrier
        Q_barrier = 1e7*eye(DynOpt.StateDim+DynOpt.nparams);
        state_diff = (b(1:DynOpt.StateDim+DynOpt.nparams)-DynOpt.OptXstory(1:DynOpt.StateDim+DynOpt.nparams,DynOpt.BackTimeIndex+i));
        
        % state bounds barrier
        % barrier on the duty cycle
        d_bound = abs(b(12));
        if b(12) > 0.1 && b(12) < 0.95
            D_barrier = 0;
        else
            D_barrier = 1e4;
        end
        
        % barrier on the amplitude
        Au_bound = abs(b(11));
        upper_bound = abs(0.99*pi/2-max(DynOpt.target_attitude));
        if b(11) > 0.05 && b(11) < upper_bound
            Au_barrier = 0;
        else
            Au_barrier = 1e4;
        end
        
        Jtot = Jtot -min(DynOpt.min_eig_past,min_eig) + D_barrier*d_bound + Au_barrier*Au_bound; %+state_diff'*Q_barrier*state_diff
        DynOpt.min_eig_past = min_eig;
    end
    
end