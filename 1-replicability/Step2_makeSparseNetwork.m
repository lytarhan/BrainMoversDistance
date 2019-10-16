% Leyla Tarhan & Evan Fields
% https://github.com/lytarhan
% 10/2019
% MATLAB R2017b

% Step 2 in use case #1 (replicability analysis): make a sparse network to
% efficiently specify the paths between voxels, along which you'll consider
% moving activation "dirt" to make one distribution of fMRI responses equal
% the other (the core computation behind the Wasserstein Distance). 

% the alternative is a fully-connected network, which takes longer to
% compute because the optimizer has to consider more possible paths. 

% input: coordinates for every voxel (i, j, k)
% output: a list of direct connections between downsampled meta-voxels. 

% To use this script: download the dataset "1-Replicability" from this project's OSF
% repository (https://osf.io/m9ac3/) to a directory called "Data", which 
% should be saved to the same directory as this script.

% To Do:
% [] add in timing estimates (for each section and at the head of the
% script).

%% clean up

clear all
close all
clc

%% File structure

% directory with distance matrix between voxels, voxel coordinates
dataDir = 'Data';

% developing over-ride:
devFlag = 1;
if devFlag
    dataDir = 'C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Project - BrainMoversDistance\Outputs\OSF - DataForGitHub\1-Replicability';
end

% directory to save the sparse edges to:
saveDir = dataDir;

addpath('../utils')

%% Setup

origVoxSize = 3.0; % how large are voxels in each dimension? (assumes isotropic)
localMaxDist = 6.0; % max distance in mm for connecting "local" voxels in step 1
% default: distance between adjacent meta-voxels (metavoxels are 2x2x2
% original voxels)

allowed_distRatio = 1.05; % initial parameter to decide when to add edges between far-apart voxels
allowed_inflation = 1.01; % parameter to decide when to delete edges 
% -- if 3 voxels are all directly connected, and the sum of 2 edges <= this factor * 3rd edge, delete the 3rd edge

distRatioStability = 0.01; % parameter to decide whether the max distance 
% ratio is "stable" between iterations of refining the network -- stop 
% when the current max distance ratio - the last iteration's max distance 
% ratio is <= this value

maxDistRatio = 1.1; % max distance ratio (captures the differences between 
% the sparse and fully-connected networks). Use this parameter to decide
% when the sparse network is done.


% load in the voxel coordinates
load(fullfile(dataDir, 'VoxelCoordinates.mat')) % coords
% coords: (i, j, k) coordinates within the Talairach brain box for every
% voxel in your subset (in the demo data, this subset is all reliable
% voxels, mainly covering high-level visual cortex)

%% meta-voxelize

% timing: a few minutes (on a standard laptop)

% downsample the data into "meta-voxels", which you'll also do with your
% data in Step 3.

[mv_to_v_mat, mv_distmat] = makeMetaVoxels(origVoxSize, coords);
disp('made meta-voxels.')

% save the meta-voxels, so you can easily load them back in during Step 3:
save(fullfile(saveDir, 'MetaVoxels.mat'), 'mv_to_v_mat', 'mv_distmat', 'coords', 'origVoxSize');
disp('saved meta-voxels.')

%% step 1: initialize the network

% timing: a few seconds

% (A) initialize a totally unconnected graph (1 vertex per meta-voxel, no edges
% connecting them):
nMetaVox = size(mv_distmat, 1);
sparseNet = sparse(nMetaVox, nMetaVox);
% store the network as a square, symmetrical sparse matrix:
    % ...0: meta-voxels aren't connected
    % ...1: meta-voxels are connected
weightedSparseNet = sparse(nMetaVox, nMetaVox);
% also set up a weighted version of this network, where each cell contains 
% an edge's distance (the straight-line distance in the brain (in mm) 
% between the voxels connected by the edge)

    
% (B) add in local edges:
% get voxels that are relatively close to each other
localVox = logical(mv_distmat <= localMaxDist); 
% get rid of the diagonal (don't connect voxels to themselves)
localVox = logical(localVox - eye(nMetaVox));
assert(all(diag(localVox)) == 0, 'Check that no voxels are connected to themselves.')
% update the sparse network (add in the local edges)
sparseNet(localVox) = 1;
% updated the weighted sparse network (add in the distances for the local
% edges)
weightedSparseNet(localVox) = mv_distmat(localVox);

fprintf('\n\nSparse Network initialized...\n')

%% step 2: add edges

% [] start again here!

% timing: c. 15-30 mins

fprintf('\nadding non-local edges...\n')

% (A) calculate the graph-distance between each pair of meta-voxels
% set up a graph
G = graph(weightedSparseNet);
% get graph distance between every pair of voxels (can only walk along
% existing edges)
graphDists = distances(G); 


% (B) add edges between voxels with a high "distance ratio"

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
distRatios = graphDists ./ mv_distmat;

% loop through the voxel pairs with high distance ratios:
for i = 1:nMetaVox
    if mod(i, 100) == 0
        fprintf('working through voxel %d / %d...\n', i, nMetaVox)
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
          sparseNet(i, j) = 1;
          sparseNet(j, i) = 1;
       end
       
    end
    
end


%% step 3: trim edges

% ...look at every "triangle" of meta-voxels (set of 3 voxels that are all
% directly connected)
    % ...go through the upper triangle of the network matrix (all i's are <
    % j)
        % ...multiply rows i & j 
        % ...look for 1's in the product (i & j were connected to k)
            % ...only look in indices >= j+1 to avoid finding the same
            % triangle multiple times
    % ...for each edge, if sum of the other 2 edges ~ that edge, delete
    % that edge 
            % ...edge1 + edge2 < allowed_inflation*edge 3 --> delete edge 3
    
% [] update distance ratios again
    
%% step 4: refine

% [] repeat steps 2 & 3 until the max distance ratio becomes stable
% [] check max distance ratio (is it within the bounds we defined above?)

% switch to adding edges more than trimming edges
% [] increase probability of adding an edge -- few ways to do this:
    % ...run adding step a few times in a row without the trimming step
    % ...increase the probabilities (multiply them all by some factor)
    % ...decrease cutting probabilities somehow?
% [] check max distance ratio again
    
% [] go through steps 2 & 3 a few more times, increasing the probability 
% of adding an edge after each iteration
% ...stop when the sparse network's efficiency is reasonably similar to the
% fully-connected network

%% save the sparse network

% format: matrix, where each row = 1 edge (col 1 = voxel 1; col 2 = voxel
% 2)

% [] save as a .csv or .mat file


    




