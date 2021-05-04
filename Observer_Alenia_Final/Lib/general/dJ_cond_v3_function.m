%% signal information policy
function [DynOpt, params] = dJ_cond_v3_function(DynOpt,params,theta,beta,gamma)
    
    if DynOpt.ActualTimeIndex > 2
    
        % position
        pos_end = DynOpt.ActualTimeIndex;
        
        % different pos_start
%         pos_start = max(1,pos_end-DynOpt.w*DynOpt.Nts);
        pos_start = 1;
        
        % temp buffers
        Y_temp = DynOpt.Y_full_story(:,pos_start:pos_end-1);
        dY_temp = DynOpt.dY_full_story(:,pos_start:pos_end-1);
        Yint_temp = DynOpt.intY_full_story(:,pos_start:pos_end-1);

        Yhat_temp = DynOpt.Yhat_full_story(:,pos_start:pos_end-1);
        dYhat_temp = DynOpt.dYhat_full_story(:,pos_start:pos_end-1);
        Yint_hat_temp = DynOpt.intYhat_full_story(:,pos_start:pos_end-1);

        % condition with error
        w = DynOpt.Nts;
        mode = 'back';
        pos = size(Y_temp,2);
        DynOpt.e = abs(sum(moving_average_single(Y_temp-Yhat_temp,w,pos,mode)));
        DynOpt.e_dot = abs(sum(moving_average_single(dY_temp-dYhat_temp,w,pos,mode)));

        initial_pos = min(w,size(Yint_temp,2)-1);
        temp_e_int = zeros(1,DynOpt.dim_out);
        for i=1:DynOpt.dim_out
            int_arg = abs((Yint_temp(i,end-initial_pos:end)-Yint_hat_temp(i,end-initial_pos:end)));
            temp_e_int(i) = trapz(DynOpt.Ts,int_arg);
        end
        DynOpt.e_int = abs(mean(temp_e_int));

        DynOpt.e_mesh = theta*DynOpt.e + beta*DynOpt.e_dot + gamma*DynOpt.e_int;
        DynOpt.dJ_cond = mean(abs(DynOpt.e_mesh));
        DynOpt.dJ_cond_story(:,end+1) = [DynOpt.e; DynOpt.e_dot; DynOpt.e_int; DynOpt.dJ_cond];
    else
        DynOpt.dJ_cond = 0;
        DynOpt.dJ_cond_story(:,end+1) = [0; 0; 0; DynOpt.dJ_cond];
    end
end