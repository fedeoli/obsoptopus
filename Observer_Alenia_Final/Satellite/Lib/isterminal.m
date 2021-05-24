%% function to get if a state is terminal
function isTerminal = isterminal(s)
    if s(1) < 1e-3 && s(2) < 1e-3
       isTerminal = 1; 
    else
       isTerminal = 0; 
    end
end