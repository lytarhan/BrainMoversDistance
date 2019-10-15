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

% directory to save the sparse edges to:
saveDir = dataDir;

addpath('../utils')

%% Setup

origVoxSize = 3.0; % how large are voxels in each dimension? (assumes isotropic)
localMaxDist = 3.0; % max distance for connecting "local" voxels in step 1
% default: distance between voxels


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

% [] load in the voxels x voxels distance matrix (?? maybe just get this from coords?)

%% meta-voxelize

% downsample the data into "meta-voxels", which you'll also do with your
% data in Step 3.
[mv_to_v_mat, mv_distmat] = makeMetaVoxels(origVoxSize, coords);
disp('made meta-voxels.')


%% step 1: initialize the network

% ...initialize a totally unconnected graph (# vertices = # meta-voxels, but 0
% edges)
% [] sparse matrix (meta-voxels x meta-voxels):
    % ...1 = voxels are connected
    % ...0 = voxels aren't connected

% ...add in local connections
    % ...connect every voxel with its nearest neighbors (defined by some
    % max distance)
    % [] get the cells in the n x n distance matrix <= localMaxDist, update the corresponding entries of the sparse network matrix 
    

%% step 2: add edges

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


    




