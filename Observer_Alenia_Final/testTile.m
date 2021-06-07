%% test tile coding
clear all
close all
clc

% bounds = [  -1 1; ...
%             -1 1; ...
%             -1 1; ...
%             -1 1; ...
%             -1 1; ...
%             -1 1; ];
% state = 0.5 - rand(1,6);
        
bounds = [  -1 1; ...
            -1 1; ];
state = [0.24, 0.34];

M = 5;
N = 4;
A = 1;

tile = build_tiles(M,N,A,bounds);


action = 1;
[feat,ind] = getFeatures(tile,M,N,A,action,state);

plot_tiles2D(tile,action,feat,ind,state);