function [lambda_min, Y_bar] = performance_index_v2(Y_story,n_rows)

    global DynOpt
    
    Y_bar = zeros(n_rows,DynOpt.Nts);
    
    for i = 1:n_rows
        for j = 1:DynOpt.Nts
           Y_bar(n_rows-i+1,DynOpt.Nts-j+1) = Y_story(1,end-(i+j)+2); 
        end
    end
    
    lambda_min = min(eig(pinv(Y_bar'*Y_bar)));

end