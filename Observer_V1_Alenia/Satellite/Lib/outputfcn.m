function stop = outputfcn(x, optimValues, state)
    global DynOpt
 
    if DynOpt.input_tuning == 0
        state_dist = abs(min(x-DynOpt.temp_x0));
    else
        state_dist = abs(min(x-DynOpt.temp_x0(end-2:end-1)));
    end
    init = (strcmp('init',state)) || (strcmp('iter',state));
    
%     if ((optimValues.fval < DynOpt.TolExit_J) || (state_dist < DynOpt.TolExit_X)) && (init==0)
    if (optimValues.fval < DynOpt.TolExit_J(1)) && (DynOpt.input_tuning == 0)
        stop = 1;
    else
        stop = 0;
    end
end