function stop = outputfcn(x, optimValues, state)
    global DynOpt
    
%     state_dist = abs(min(x-DynOpt.temp_x0));
    init = strcmp('init',state);    


%     if ((optimValues.fval < DynOpt.TolExit_J) || (state_dist < DynOpt.TolExit_X)) && (init==0)
    if optimValues.fval < DynOpt.TolExit_J(1) &&  ~init
        stop = 1;
        DynOpt.J_big = 0;
    elseif (optimValues.fval > DynOpt.TolExit_J(2)) && init  
        stop = 0;
        DynOpt.J_big = 1;
    else
        DynOpt.J_big = 0;
        stop = 0;
    end
    
    % uncomment to disable the function
    stop = 0;
end