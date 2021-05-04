%% function plot Q
function Qval = plotQ(Q)
    %%% show results
    figure(1)
    clf
    a = max(Q,[],3);
    Qval = -1*ones(size(a,1),size(a,2));
    for j=1:size(a,1)
        for k=1:size(a,2)
           b = find(Q == a(j,k));
           [~,~,Aval] = ind2sub(size(Q),b);
           Qval(j,k) = Aval;
        end
    end
    bar3(Qval)
    xlabel('action')
    ylabel('state')
    zlabel('value')
end