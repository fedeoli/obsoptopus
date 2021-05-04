function stop = outputfcn(x, optimValues, state)
 
    if (optimValues.fval < 1e-8)
        stop = 1;
    else
        stop = 0;
    end
end