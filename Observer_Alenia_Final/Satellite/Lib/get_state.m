%% get state from DynOpt 
function nu = get_state(DynOpt)

    % get data
    edot = DynOpt.dJ_cond_story(2,:);
    eint = DynOpt.dJ_cond_story(3,:);
    
    % L2 norm
    L2_edot = trapz(DynOpt.Ts,edot.^2);
    L2_eint = trapz(DynOpt.Ts,eint.^2);
    
    % get state
    nu = [L2_edot, L2_eint];
end