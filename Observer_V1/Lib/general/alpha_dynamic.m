function alpha = alpha_dynamic(state)

    global DynOpt
    
    if DynOpt.alpha_dyn == 1
        alpha = DynOpt.alpha_grad.*abs(mean(state));
    else
        alpha = DynOpt.alpha_grad;
    end
   
end