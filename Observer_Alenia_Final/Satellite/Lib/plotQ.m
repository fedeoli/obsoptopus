%% function plot Q
function Qval = plotQ(Q,domain_S)

    
    %%% show results
    figure(1)
    clf
    a = max(Q,[],2);
    Qval = -1*ones(size(a,1),1);
    for j=1:size(a,1)
       Aval = find(Q(j,:) == a(j));
       Qval(j) = Aval;
    end
    
    size_Q = [size(domain_S,2),size(domain_S,2)];
    Qval = reshape(Qval,size_Q);
%     bar3(Qval)
    heatmap(Qval);
    xlabel('e')
    ylabel('edot')
%     zlabel('value')
end