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

% input: voxels x voxels matrix containing the real distances between every
% pair of voxels (e.g., in mm). 
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

nLocalEdges = 6; % how many local connections to make in the first step
origVoxSize = 3.0; % how large are voxels in each dimension? (assumes isotropic)

pAddEdges = 1; % [] 
% specify the probability for adding edges during step 2
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

% ...initialize a totally unconnected graph (# vertices = # voxels, but 0
% edges)

% ...add in local connections
    % ...connect every voxel with its 6 nearest neighbors
    

%% step 2: add edges

% ...calculate the distance ratio for each pair of voxels
    % ...graph-distance / actual brain distance
    % ...1 edge in graph-space = length of 1 (?)
    % ...if the value of the distance ratio is much greater than 1, the
    % sparse network is much less efficient than the fully-connected one
    % (you have to take a circuitous path to get between nearby voxels)
    
% ...add an edge between each pair of meta-voxels according to some
% probability function
% [] specify the probability function
% ...voxel pairs with higher distance ratios are more likely to be
% connected by a new edge
% [] update the distance ratios for all the voxel pairs

%% step 3: trim edges

% ...look at every "triangle" of meta-voxels (set of 3 voxels that are all
% directly connected)
    % ...for each edge, if sum of the other 2 edges ~ that edge, delete
    % that edge
    
%% step 4: refine

% [] repeat steps 2 & 3 until the max distance ratio becomes stable
% [] increase probability of adding an edge

% [] go through steps 2 & 3 a few more times, increasing the probability 
% of adding an edge after each iteration
% ...stop when the sparse network is reasonably similar to the
% fully-connected network


    




