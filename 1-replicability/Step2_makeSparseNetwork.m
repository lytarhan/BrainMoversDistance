% Leyla Tarhan & Evan Fields
% https://github.com/lytarhan
% 10/2019
% MATLAB R2017b

% *** THIS CODE IS STILL IN PROGRESS ***

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

clc
fprintf('Step 0: Downsampling the Data\n')
fprintf('...................................\n\n')

% downsample the data into "meta-voxels", which you'll also do with your
% data in Step 3.
[mv_to_v_mat, mv_distmat] = makeMetaVoxels(origVoxSize, coords);
disp('made meta-voxels.')

% save the meta-voxels, so you can easily load them back in during Step 3:
save(fullfile(saveDir, 'MetaVoxels.mat'), 'mv_to_v_mat', 'mv_distmat', 'coords', 'origVoxSize');
disp('saved meta-voxels.')

%% step 1: initialize the network

% timing: a few seconds

clc
fprintf('Step 1: Initializing the Network\n')
fprintf('...................................\n\n')

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
% update the weighted sparse network (add in the distances for the local
% edges)
weightedSparseNet(localVox) = mv_distmat(localVox);

fprintf('\n\nSparse Network initialized...\n')

%% step 2: add edges

% timing: c. 1 min.


clc
fprintf('Step 2: Adding Edges\n')
fprintf('........................\n\n')

sparseNet2 = addEdges(sparseNet, weightedSparseNet, mv_distmat, allowed_distRatio);


%% step 3: trim edges

% timing: % [] 

clc
fprintf('Step 3: Trimming Edges\n')
fprintf('........................\n\n')

[sparseNet3, weightedSparseNet2] = trimEdges(sparseNet2, mv_distmat, allowed_inflation);

% get the new max distance ratio:
currMaxRatio = getMaxDistRatio(weightedSparseNet2, mv_distmat);
    
%% step 4: stabilize the network
% repeat adding and trimming steps until the network becomes stable

% timing: % [] 


clc
fprintf('Step 4: Stabilizing the Network\n')
fprintf('................................\n\n')


distRatioDiff = 1; % initialize the difference between max distance ratios 
% across iterations to something arbitrary but high
prevMaxRatio = currMaxRatio;
sparseNet4b = sparseNet3;
weightedSparseNet4b = weightedSparseNet2;

while distRatioDiff > distRatioStability
   % add edges:
   sparseNet4a = addEdges(sparseNet4b, weightedSparseNet4b, mv_distmat, allowed_distRatio);
    
   % trim edges:
   [sparseNet4b, weightedSparseNet4b] = trimEdges(sparseNet4a, mv_distmat, allowed_inflation);
   
   % re-evaluate max distance ratio:
   currMaxRatio = getMaxDistRatio(weightedSparseNet4b, mv_distmat);
   
   % get difference from the last iteration
   distRatioDiff = prevMaxRatio - currMaxRatio;
   prevMaxRatio = currMaxRatio; % update for the next iteration
end

%% step 5: prioritize adding edges

% timing: % [] 


clc
fprintf('Step 5: Refining the Network\n')
fprintf('................................\n\n')


% check max distance ratio (is it within the bounds we defined above?)
if currMaxRatio <= maxDistRatio
    disp('Sparse Network is DONE!')
else % switch to adding edges more than trimming edges
    fudgeFactor = 0.01; % allowable difference from ideal distance ratio
    currAllowed_distRatio = allowed_distRatio; % parameter to determine probability of adding edges
    sparseNet5c = sparseNet4b;
    weightedSparseNet5c = weightedSparseNet4b;
    
    while maxDistRatio < currMaxRatio - fudgeFactor
       % increase the probability of adding edges
       currAllowed_distRatio = currAllowed_distRatio - 0.1; % smaller value = more likely to add an edge
       
       % add edges:
       fprintf('\nadding edges - step 1:\n')
       sparseNet5a = addEdges(sparseNet5c, weightedSparseNet5c, mv_distmat, currAllowed_distRatio);
       
       % add edges again:
       fprintf('\nadding edges - step 2:\n')
       sparseNet5b = addEdges(sparseNet5a, weightedSparseNet5c, mv_distmat, currAllowed_distRatio);
       
       % trim edges:
       [sparseNet5c, weightedSparseNet5c] = trimEdges(sparseNet5b, mv_distmat, allowed_inflation);
       
       % re-evaluate max distance ratio:
       currMaxRatio = getMaxDistRatio(weightedSparseNet5c, mv_distmat);
    
       
    end
    disp('Sparse Network is DONE!')
end


%% save the sparse network

% format: matrix, where each row = 1 edge (col 1 = voxel 1; col 2 = voxel
% 2)
[Vox1, Vox2] = find(sparseNet5c); % get indices for the non-zero elements in the sparse net (symmetrical)
T = table(Vox1, Vox2);

% write to .csv:
disp('Saving sparse network...')
writetable(T, fullfile(saveDir, 'sparseEdges.csv'), 'Delimiter', ',');
fprintf('\n\nDONE!! \nSaved sparse network as a .csvin dataDir.\n')


    




