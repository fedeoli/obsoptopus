function alpha = alpha_dynamic(state)

    global DynOpt
    
    alpha = 1e-6.*abs(state);
%     alpha = DynOpt.alpha_grad;

    
end