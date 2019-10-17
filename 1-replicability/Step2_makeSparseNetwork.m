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

% timing: c. 1 min.

sparseNet2 = addEdges(sparseNet, weightedSparseNet, mv_distmat, allowed_distRatio);


%% step 3: trim edges

% timing: % [] 

% [] start again here (takes awhile)
sparseNet3 = trimEdges(sparseNet2, mv_distmat, allowed_inflation);

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


    




