function stop = outputfcn(x, optimValues, state)
    global DynOpt
 
    state_dist = abs(min(x-DynOpt.temp_x0));
    init = (strcmp('init',state)) || (strcmp('iter',state));
    
%     if ((optimValues.fval < DynOpt.TolExit_J) || (state_dist < DynOpt.TolExit_X)) && (init==0)
    if (optimValues.fval < DynOpt.TolExit_J(1))
        stop = 1;
    else
        stop = 0;
    end
end