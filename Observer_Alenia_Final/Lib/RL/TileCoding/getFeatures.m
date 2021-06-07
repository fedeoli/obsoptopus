%% function to localise feature
function [feat,ind] = getFeatures(tile,M,N,A,action,state)
    
    % init indices
    ndim = length(state);
    ind = 0;    
    feat = zeros(A*N*(M+1)^ndim,1);
    
    % cycle over the actions
    for a = 1:A
        % cycle over tiles
        for n=1:N

           if a==action
               % get grid
               grid = tile(action,n).grid;

               % check grid position
               sub = zeros(ndim,1);
               str = [];
               str_state = [];
               for d = 1:ndim
                   temp = find(state(d) >= grid(d,:),1,'last');
                   if ~isempty(temp)
                       sub(d) = min(temp,M+1);
                   else
                       sub(d) = 0;
                   end
                   sub_str = sprintf('sub(%d),',d);
                   str = strcat(str,sub_str);
                   str_state = strcat(str_state,sprintf('M+1,'));
               end
               
               if nnz(sub) == ndim
                    try
                    eval(strcat('ind(end+1) = sub2ind([A,N,',str_state(1:end-1),'],a,n,',str(1:end-1),');'));
                    catch
                       fail = 1; 
                    end
               end
           end
        end
    end
    
    % assign feature
    ind = ind(2:end);
    feat(ind) = 1;

end