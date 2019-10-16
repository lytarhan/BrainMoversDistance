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

% downsample the data into "meta-voxels", which you'll also do with your
% data in Step 3.

% N.B.: this step takes a few minutes on a standard laptop
[mv_to_v_mat, mv_distmat] = makeMetaVoxels(origVoxSize, coords);
disp('made meta-voxels.')

% save the meta-voxels, so you can easily load them back in during Step 3:
save(fullfile(saveDir, 'MetaVoxels.mat'), 'mv_to_v_mat', 'mv_distmat', 'coords', 'origVoxSize');
disp('saved meta-voxels.')

%% step 1: initialize the network

% (A) initialize a totally unconnected graph (1 vertex per meta-voxel, no edges
% connecting them):
nMetaVox = size(mv_distmat, 1);
sparseNet = sparse(nMetaVox, nMetaVox);
% store the network as a square, symmetrical sparse matrix:
    % ...0: meta-voxels aren't connected
    % ...1: meta-voxels are connected
weightedSparseNet = sparse(nMetaVox, nMetaVox);
% also set up a weighted sparse net (each cell = distance that the edge
% traverses between the 2 meta-voxels)

    
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


%% step 2: add edges

% (A) calculate the graph-distance between each pair of meta-voxels
% set up a graph
G = graph(weightedSparseNet);

% [] start again here! (weights are in G.Edges table already)
% d = distances(G) returns a matrix, d, where d(i,j) is the length of the 
% shortest path between node i and node j. If the graph is weighted 
% (that is, G.Edges contains a variable Weight), then those weights are 
% used as the distances along the edges in the graph. Otherwise, all edge 
% distances are taken to be 1.


% ...calculate the distance ratio for each pair of voxels
    % ...graph-distance / actual brain distance
    % ...length of each edge = actual brain distance between the connected voxels
    % ...if the value of the distance ratio is much greater than 1, the
    % sparse network is much less efficient than the fully-connected one
    % (you have to take a circuitous path to get between nearby voxels)
    % [] find a function like floyd_warshall_shortest_paths that find the
    % shortest path between voxels
    
% ...add an edge between each pair of meta-voxels according to some
% probability function (IF they're not already directly connected)
        % ...if distance ratio <= 1.05, probability of adding
        % an edge = 0 (don't add an edge)
        % ...if distance ratio > 1.05, probability of adding an edge =
        % distance ratio - 1.05 (could be > 1)
        % ...decide whether to add this edge:
            % ...generate a random # [0, 1]. If that random # is < the probability of adding an edge, add the edge. Otherwise, don't add it. 
        
% ...voxel pairs with higher distance ratios are more likely to be
% connected by a new edge


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


    




