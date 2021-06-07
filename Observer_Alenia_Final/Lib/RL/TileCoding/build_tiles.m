%% build tiles for tile coding (RL)
function tile = build_tiles(M,N,A,bounds)
    
    % unpack bounds
    ndim = size(bounds,1);
    
    % create grids
    grid = zeros(ndim,M+1);
    for i=1:ndim
        grid(i,:) = linspace(bounds(i,1), bounds(i,2),M+1);
    end
    
    % create tiles
    tile = [];
    width = diff(grid(:,1:2),1,2);
    offset = width/N;
    shift = 1:2:2*ndim;
    
    for a = 1:A
       for n = 1:N
            newgrid = zeros(ndim,M+1);
            for d = 1:ndim
               newgrid(d,:) = grid(d,:) + (n-1)*shift(d)*offset(d); 
            end
            tile(a,n).grid = newgrid;
       end
    end
end