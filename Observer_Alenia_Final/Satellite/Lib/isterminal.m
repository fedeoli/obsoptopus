%% function to get if a state is terminal
function RL = isterminal(s,terminal_cond,terminal_streak_cond,RL)
    if (sum(abs(s) < terminal_cond) == length(s))
        RL.S.terminal_streak = RL.S.terminal_streak + 1;
    else
        RL.S.terminal_streak = 0;
    end
    
    if RL.S.terminal_streak > terminal_streak_cond
        RL.S.isTerminal = 1;
    else
        RL.S.isTerminal = 0;
    end
end