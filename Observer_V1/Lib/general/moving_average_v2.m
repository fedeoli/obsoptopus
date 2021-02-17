function y_out = moving_average_v2(y_in,w,mode)

    % init 
    n_dati = length(y_in); 
    y_out = zeros(n_dati,1);
    
    
    for i=1:n_dati
        pos = i;
        % moving average
        if (strcmp(mode,'center') == 1)
            if (pos-w > 0) && ((pos+w) <= n_dati)
                y_out(pos,1) = mean(y_in(pos-w:pos+w,:));
            else
                y_out(pos,1) = y_in(pos,:);
            end
        elseif (strcmp(mode,'back') == 1)
            if (pos-w > 0)
                y_out(pos,1) = mean(y_in(:,pos-w:pos),2);
            end
        elseif (strcmp(mode,'forward') == 1)
            if (pos-w > 0)
                y_out(pos,1) = mean(y_in(:,pos:pos+w),2);
            end
        end
        
    end

end