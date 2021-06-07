%% function to get if a state is terminal
function isTerminal = isterminal(s,terminal_cond)
    if (sum(abs(s) < terminal_cond) == length(s))
        isTerminal = 1;
    else
        isTerminal = 0; 
    end
end