%% test tile coding
function tile = init_tile(statedim,A,M,N)

    bounds = [0 1].*ones(statedim,2);
    tile = build_tiles(M,N,A,bounds);

end