function set_input(pos)

    global DynOpt params
    
    
    % system control input
    if DynOpt.params.input_flag == 1
        u = params.U(pos);
        params.u = u;
        DynOpt.params.u = u;
    else
        params.u = 0;
        DynOpt.params.u = 0;
    end
    
end