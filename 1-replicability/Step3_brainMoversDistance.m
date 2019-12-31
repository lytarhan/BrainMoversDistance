% Leyla Tarhan
% https://github.com/lytarhan
% 10/2019
% MATLAB R2017b


% Step 3 in use case #1 (replicability analysis): calculate the wasserstein
% distance between all of the splits of the data, for each condition.

% on each iteration, calculate the wasserstein distance between two
% "brains": group 1 (averaged neural responses across 6/13 subjects) and
% group 2 (averaged neural responses across 7/13 subjects). Do this
% separately for each condition in the data (in this case, 1 condition = 1
% video clip).

% at this point, the data have already been restricted to *reliable* voxels
% using the method demonstrated here: https://osf.io/m9ykh/.

% to get an estimate of the baseline distance between dtaa splits, set
% "scramFlag" = 1 (randomly permutes the order of the conditions in one
% group, which should result in a very high distance between groups). 

% -------------------------------------------------------------------------
% To use this script:
% (1) install the Gurobi optimizer on your machine 
% (https://www.gurobi.com/documentation/quickstart.html). An academic 
% license can be obtained at no cost. 

% (2) download the dataset "1-Replicability" from this project's OSF repository 
% (https://osf.io/m9ac3/) to a directory called "Data", which 
% should be saved to the same directory as this script.


% NB: This script takes a long time to run (up to 24 hours on my PC). 
% I suggest letting it run over the weekend.


%% Clean up

clear all
close all
clc

%% file structure

dataDir = fullfile('Data');
saveDir = 'Results';
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% developing over-ride:
devFlag = 1;
if devFlag
    dataDir = 'C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Project - BrainMoversDistance\Outputs\OSF - DataForGitHub\1-Replicability';
end


addpath('../utils')

%% set up some parameters

% scrambling to get a baseline?
scramFlag = 0;

% ran step 2? (default = 1):
step2Flag = 1;

% how many splits in the data? (must match # of splits made in
% Step1_makeSplits)
nSplits = 100;

% how many conditions in your data (1 condition = 1 video clip)?
nConds = 120;

% size of the original voxels (assuming they're isotropic)
origVoxSize = 3.0; % mm

if scramFlag
    saveDir = 'Results-Scrambled';
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    
    disp('starting analysis on scrambled data...')
else
    disp('starting analysis on intact data...')
end

%% Load in general data

% voxel coordinates
load(fullfile(dataDir, 'VoxelCoordinates.mat')) % coords

% sparse network specifying a subset of edges between meta-voxels (made in step 2):
sn = readtable(fullfile(dataDir, 'sparseEdges.csv'));

whos

%% Make meta-voxels

% either downsample the data here or load the meta-voxels made in step 2:
if step2Flag
    load(fullfile(dataDir, 'MetaVoxels.mat'));
else
    % downsample the data into "meta-voxels" (2x2x2 original voxels).
    [mv_to_v_mat, mv_distmat] = makeMetaVoxels(origVoxSize, coords);
end
disp('made meta-voxels.')

% output: 
% (1) sparse matrix (original voxels x meta-voxels) marking members of each
% meta-voxel
% (2) distance matrix for the meta-voxels

% set up for the model:
[A, cTranspose] = makeSparseConstraintMatrix(mv_distmat, sn);
disp('made sparse A and c-transpose')

% if you don't have a sparse network, just do this for all edges (voxel pairs):
% [A, cTranspose] = makeConstraintMatrix(mv_distmat);
% disp('made A and c-transpose')

%% Loop through the splits

tic
for s = 1:nSplits
    %% Load in the data

    % activations for current split: group 1 and group 2
    fn = sprintf('groupAverages_split%d.mat', s);
    load(fullfile(dataDir, 'split averages', fn)) % group1Avg & group2Avg
    assert(size(group1Avg, 2) == nConds, 'try transposing data matrices.')
    
    % option to scramble:
    if scramFlag
       nVox = size(group2Avg, 1); 
       scramIdx = randperm(nVox);
       group2Avg = group2Avg(scramIdx, :);        
    end
        
    %% calculate wasserstein distance for each condition

    % set up vector to store distances (1/condition):
    wd = zeros(nConds, 1);
    for c = 1:nConds % loop through conditions
        clc
        fprintf('split %d / %d: calculating wasserstein distance for condition %d / %d...\n', s, nSplits, c, nConds);
        
        % pull out data for current condition:
        group1Curr = group1Avg(:, c); % nVoxels x 1
        group2Curr = group2Avg(:, c);
        
        %% format the data
        
        disp('formatting the data:')
        
        % normalize:
        [g1Norm, g2Norm] = normalize_wasserstein(group1Curr, group2Curr);
        % output: 2 vectors (normalized activations for group 1 and group 2)
        disp('...normalized')
        
        % downsample to "meta-voxels":
        [g1MV, g2MV] = metavoxelize(g1Norm, g2Norm, mv_to_v_mat);
        % output: 2 vectors (downsampled activations for group 1 and group 2)
        disp('...down-sampled.')
        
        
        %% wasserstein distance
        wd(c) = computeWasserstein(g1MV, g2MV, A, cTranspose);
        % output: avg. cost (in mm of brain) to move one unit of normalized
        % activation from group 1 to group 2 for this condition.
        
        
    end
    
    
    %% save the output
    % save the wasserstein distances for this split
    dn = sprintf('wasserstein_split%d.mat', s);
    save(fullfile(saveDir, dn), 'wd') % 1 wasserstein distance per condition
end

endTime = toc;
fprintf('calculated wasserstein distances for all conditions across %d splits in %0.2f hours!\n', nSplits, endTime/3600)






