function y_out = moving_average(y_in,w,pos,mode)

    % init 
    y_out = y_in;
    
    % moving average
    if (strcmp(mode,'center') == 1)
        if (pos-w > 0) && ((pos+w) <= length(y_in))
            y_out(:,pos) = mean(y_in(:,pos-w:pos+w),2);
        end
    elseif (strcmp(mode,'back') == 1)
        if (pos-w > 0)
            y_out(:,pos) = mean(y_in(:,pos-w:pos),2);
        end
    elseif (strcmp(mode,'forward') == 1)
        if (pos-w > 0)
            y_out(:,pos) = mean(y_in(:,pos:pos+w),2);
        end
    end

end