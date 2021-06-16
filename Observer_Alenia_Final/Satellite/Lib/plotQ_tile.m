%% 
function plotQ_tile(RL)

    bounds = [0 1].*ones(RL.S.statedim,2);
    [W,Q] = meshgrid(bounds(1,1):0.01:bounds(1,2),bounds(2,1):0.01:bounds(2,2));
    NW = prod(size(W));
    NQ = prod(size(Q));
    
    Q_handle = @(w,q)Qfun(w,q,RL);
    
    Q = Q_handle(W,Q);
    
    mesh(NW,NQ,Q)
    
end