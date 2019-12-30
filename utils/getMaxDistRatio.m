function currMaxRatio = getMaxDistRatio(weightedNet, voxDists)
% Leyla Tarhan & Evan Fields
% get the maximum distance ratio between the sparse network and the
% fully-connected one. This gives you an idea of how much *less* efficient
% the sparse network is than the fully-connected one.

% inputs:
    % ...weightedNet: nVoxels x nVoxels. If 2 voxels are directly connected 
    % by an edge, the corresponding cell = actual distance (in the brain)
    % between those 2 voxels.
    % ...voxDists: nVoxels x nVoxels. Each cell = actual distance (in the
    % brain, mm) between voxels i & j.
    
% output:
    % ...currMaxRatio: the maximum distance ratio of all the edges.
    
%% setup

fprintf('\nevaluating the network''s efficiency...\n')

% set up a graph object from the weighted network
G = graph(weightedNet);

%% get distance ratio for every pair of voxels

% get graph distance between every pair of voxels (can only walk along
% existing edges)
graphDists = distances(G); 

% get distance ratio for every pair of voxels:
distRatios = graphDists ./ voxDists;

%% get max distance ratio

currMaxRatio = max(max(distRatios));


end