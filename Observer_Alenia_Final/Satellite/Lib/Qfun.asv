%%
function Q = Qfun(w,q,RL)

    for i = 1:RL.A.dimActions
        [x(:,i),x_ind(i)] = getFeatures(RL.S.tile,RL.S.M,RL.S.N,RL.A.dimActions,i,[w,q]);
        V(i) = transpose(RL.S.w)*x(:,i);
    end
    
    Q = max(V);
end