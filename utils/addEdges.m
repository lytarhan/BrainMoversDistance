function updatedNet = addEdges(sparseNet, weightedNet, voxDists, allowed_distRatio)
% Leyla Tarhan & Evan Fields
% add connections between voxels that are connected too circuitously in the
% sparse network.

% inputs:
    % ...sparseNet: nVoxels x nVoxels. If 2 voxels are directly connected
    % by an edge, the corresponding cell = 1.
    % ...weightedNet: nVoxels x nVoxels. If 2 voxels are directly connected 
    % by an edge, the corresponding cell = actual distance (in the brain)
    % between those 2 voxels.
    % ...voxDists: nVoxels x nVoxels. Each cell = actual distance (in the
    % brain, mm) between voxels i & j.
    % ...allowed_distRatio: parameter to decide when to add edges between
    % far-apart voxels. Should be >1 but not by much.
    
% output:
    % ...updatedNet: updated version of sparseNet (with more edges).


%% setup

fprintf('\nadding edges...\n')

% how many voxels?
nVox = size(voxDists, 1);

% set up a graph object from the weighted network
G = graph(weightedNet);

% set up for the new sparse network (stored as a sparse matrix -- much 
% faster to store locations and values of non-zero elements than to 
% store/update a sparse matrix in the loop below)
edgeVals = []; % 1 for every edge
edgeRow = []; % voxel i for every edge
edgeCol = []; % voxel j for every edge

%% calculate the graph-distance between each pair of meta-voxels

% get graph distance between every pair of voxels (can only walk along
% existing edges)
graphDists = distances(G); 

%% add edges between voxels with a high "distance ratio"

% distance ratio: graph-distance / actual brain distance -- this measure
% captures how much more efficient a fully-connected network is, compared
% to the sparse network at its current state. Eventually, the sparse
% network's distances should closely resemble those in full network 
% (it shouldn't have too many long, circuitous connections between voxels),
% but with far fewer edges. If the distance ratio between 2 voxels is much 
% larger than 1, the sparse network's connection between the voxels is much
% less efficient than in the fully-connected one (and we should add more 
% edges to connect these voxels directly).

% get distance ratio for every pair of voxels:
distRatios = graphDists ./ voxDists;

% loop through the voxel pairs with high distance ratios:
tic
for i = 1:nVox
    if mod(i, 100) == 0
        fprintf('working through voxel %d / %d...\n', i, nVox)
    end
    
    % get the voxels with which this voxel has a high distance ratio:
    highDistVox = find(distRatios(i, :) > allowed_distRatio);
    
    % loop through those voxels:
    for v = 1:length(highDistVox)
       j = highDistVox(v); 
       
       % pull out the distance ratio between voxels i & j
       currDistRatio = distRatios(i, j);
        
       % probabilistically add an edge between voxels i & j, so that you're
       % more likely to add an edge between voxels with higer distance
       % ratios
       prEdge = currDistRatio - allowed_distRatio;
       if rand(1) < prEdge && sparseNet(i, j) ~= 1
          % add the edge (symmetrically)
          edgeVals = [edgeVals; 1]; % edge
          edgeRow = [edgeRow; i]; % voxel i
          edgeCol = [edgeCol; j]; % voxel j
          
          % make it symmetrical:
          edgeVals = [edgeVals; 1]; % edge
          edgeRow = [edgeRow; j]; % voxel j
          edgeCol = [edgeCol; i]; % voxel i
          
       end
       
    end
    
end
fprintf('...finished adding edges in %d sec!\n', toc)

%% Make the updated sparse network

% add in the edges that already existed:
[prevRow, prevCol] = find(sparseNet); % get indices for the non-zero elements in the sparse net (symmetrical)
prevVals = ones(length(prevRow), 1); % array to store edges (1's if edge exists)

% add to the vectors we've been building up:
edgeVals = [edgeVals; prevVals];
edgeRow = [edgeRow; prevRow];
edgeCol = [edgeCol; prevCol];

% make the final new sparseNet.
updatedNet = sparse(edgeRow, edgeCol, edgeVals, nVox, nVox);

% get rid of potential duplicates:
updatedNet = spones(updatedNet);


end