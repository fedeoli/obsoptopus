%% function to locate state in RL algorithm
function [pos_final, pos, state] = locate_state(nu,domain_S)

    % state dim
    dim = size(domain_S,1);
    pos = zeros(dim,1);
    state = -1*ones(dim,1);
    
    % iterate
    for i=1:dim
        temp = find(domain_S(i,:) >= nu(i),1,'first');
        if isempty(temp)
            pos(i) = length(domain_S(i,:));
        else
            pos(i) = max(1,temp);
        end
        state(i) = domain_S(i,pos(i));
    end
    
    % from readable pos to index
    size_matrix = [size(domain_S,2),size(domain_S,2)];
    pos_final = sub2ind(size_matrix,pos(1),pos(2));
end