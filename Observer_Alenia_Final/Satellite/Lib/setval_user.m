%% function to get val from index and matrix 
function policy = setval_user(policy,pos,val)
    policy(pos(1),pos(2)) = val;
end