%% function to plot 2D tiles
function plot_tiles2D(tile,action,feat,ind,state)

    if length(state) > 2
        disp('state too big')
        return
    end
    
    figure
    hold on
    grid on
    
    lbx = Inf;
    lby = Inf;
    ubx = -Inf;
    uby = -Inf;
    
    A = size(tile,1);
    N = size(tile,2);
    M = size(tile(1,1).grid,2)-1;
    
    for i=1:N
        color(i,:) = [rand,rand,rand];
    end
    
    for i = 1:N
       
       
        
       horzval = tile(1,i).grid(1,:);
       vertval = tile(1,i).grid(2,:);
       
       lbx = min(horzval(1),lbx);
       lby = min(vertval(1),lby);
       ubx = max(horzval(end),ubx);
       uby = max(vertval(end),uby);
       
       horzlen = length(horzval);
       for j = 1:horzlen
           plot(horzval(j)*ones(1,horzlen),vertval,'--','Color',color(i,:),'LineWidth',2);
       end
       
       vertlen = length(vertval);
       for j = 1:vertlen
           plot(horzval,vertval(j)*ones(1,vertlen),'--','Color',color(i,:),'LineWidth',2);
       end
    end
    
    %%% find intersections
    nind = length(ind);    
    statedim = length(state);
    sub = zeros(N,4);
    for i = 1:nind
        
        str = '[';
        for s = 1:4
           str_temp = sprintf('sub(i,%d),',s);
           str = strcat(str,str_temp);
        end
        str = strcat(str(1:end-1),']');
    
        eval(strcat(str,' = ind2sub([A,N,M+1,M+1],ind(i));'));        
    end
    
    % get cells
    bounds = zeros(nind,2*statedim);
    for i=1:nind
        bounds(i,1:2) = tile(action,sub(i,2)).grid(1,sub(i,3):sub(i,3)+1);
        bounds(i,3:4) = tile(action,sub(i,2)).grid(2,sub(i,4):sub(i,4)+1);
    end

    for i=1:nind
        square_x = [bounds(i,1) bounds(i,2) bounds(i,2) bounds(i,1)];
        square_y = [bounds(i,3) bounds(i,3) bounds(i,4) bounds(i,4)];
        color_fill = color(i,:);
        temp = fill(square_x, square_y, color_fill,'LineWidth',2);
        temp.FaceAlpha = 0.2;
    end
    plot(state(1),state(2),'k*','LineWidth',2)
    
    xlim([lbx ubx]);
    ylim([lby uby]);
end