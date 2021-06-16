%% signal information policy
function DynOpt = dJ_cond_v5_function(DynOpt,params)
    
    buffer_ready = (nnz(DynOpt.Y_space) >= 1) || (size(DynOpt.Yhat_full_story,2) > DynOpt.Nts);

    if buffer_ready
        
        [~, pos_hat] = find(DynOpt.Y_space ~= 0);
        pos_hat = DynOpt.Y_space(pos_hat)-1;
        Yhat = DynOpt.Yhat_full_story(:,pos_hat);
        Y_buf = DynOpt.Y_full_story(:,pos_hat);
        bias_hat = zeros(9,1);
        bias_hat(1:length(params.bias)) = params.bias;
        
        n_sum = size(Y_buf,2);
        temp_e = 0;
        for i=1:n_sum
            temp_e = temp_e + norm(Y_buf(:,i)-(Yhat(:,i)+bias_hat));
        end
        
        DynOpt.e = temp_e;
        DynOpt.e_dot = 0;
        DynOpt.e_int = 0;
        
        DynOpt.dJ_cond = norm(DynOpt.e);
        DynOpt.dJ_cond_story(:,end+1) = [DynOpt.e; DynOpt.e_dot; DynOpt.e_int; DynOpt.dJ_cond];
    else
        DynOpt.dJ_cond = 1;
        DynOpt.dJ_cond_story(:,end+1) = [0; 0; 0; DynOpt.dJ_cond];
    end
end