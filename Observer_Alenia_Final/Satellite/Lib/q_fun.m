%% value function - approximation
function q = q_fun(s,a,w)

    % feature
    x = [s; a];
    
    % model
    q = transpose(w)*x;
end