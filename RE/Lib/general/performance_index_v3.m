function [lambda_min, Y_bar] = performance_index_v3(Y_story)

    global DynOpt
    
    n_rows = length(Y_story)-DynOpt.w;
    Y_bar = zeros(n_rows,DynOpt.w);
    
    
    for i = 1:n_rows
        for j = 1:DynOpt.w
           Y_bar(n_rows-i+1,DynOpt.w-j+1) = Y_story(1,end-(i+j)+2); 
        end
    end
    
    lambda_min = min(eig(pinv(Y_bar'*Y_bar)));

end