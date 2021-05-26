%% function to get if a state is terminal
function isTerminal = isterminal(s)
    terminal_cond = [5e-4*ones(3,1); 1e-4*ones(3,1)];
    if (sum(abs(s) < terminal_cond) == length(s))
        isTerminal = 1;
    else
        isTerminal = 0; 
    end
end