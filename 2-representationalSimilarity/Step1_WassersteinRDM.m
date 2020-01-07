% Leyla Tarhan
% https://github.com/lytarhan
% 1/2020
% MATLAB R2017b


% Step 1 in use case #2 (representational similarity analysis): calculate 
% the neural dissimilarities between all items in an fMRI experiment, using
% Wasserstein Distance. Do this separately for each subject, within a
% pre-defined region of interest. 

% timing notes: in an example subject with about 500 voxels, takes about 
% .21 mins to calculate 100 distances. So, in the demo data (with 72
% conditions, or 2556 pairs of conditions, it took about 8.9 hours per subject).

% N.B.: you could save time by making a sparse network first, like in use
% case #1. However, since ROI's are generally smaller than the large swathe
% of cortex that I was working with in use case #1, this isn't as necessary
% here.

% -------------------------------------------------------------------------
% To use this script:
% (1) install the Gurobi optimizer on your machine 
% (https://www.gurobi.com/documentation/quickstart.html). An academic 
% license can be obtained at no cost. 

% (2) make sure you have the parallel computing MATLAB toolbox installed.

% (3) download the dataset "2-Replicability" from this project's OSF repository 
% (https://osf.io/m9ac3/) to a directory called "Data", which 
% should be saved to the same directory as this script. This data contains 
% two parts:
    % (1) fMRI data (matrix of beta weights within a particular region for
    % each subject. The dimensions of this matrix are # voxels in the 
    % region x # stimuli (objects).
    
    % (2) group-level behavioral data from a visual search experiment,
    % involving the same stimuli. The dimensions here are # stimuli x #
    % stimuli, where each cell (i, j) in the matrix contains the 
    % log-transformed reaction time when subjects saw target i among
    % distractors j.
    
% If this dataset is not available on OSF, format your data in a similar 
% fashion. Make a struct for each subject, with the following key
% variables:

    % -betas: neural responses in the region of interest (voxels x items)
    % -voxelCoords: 3-d coordinates for every voxel within a Talairach cube (voxels x 3).
    % -voxelDists: distance matrix containing the actual distance between
    % every pair of voxels (in mm; I used straight-line distance but you
    % could also measure geodesic distance, connectivity distance, etc.) (#
    % voxels x # voxels)
    % -corrRDM: comparison neural RDM, based on another distance metric
    % (e.g., correlation) -- this isn't necessary until Step 2.
    
% save all of these subject-specific structs in a single .mat file named
% 'FormattedData_allSubs.mat' in a directory at the same level as this
% script, called 'Data'.



%% clean up
clear all
close all
clc

%% file structure

dataDir = 'Data-fMRI';
saveDir = fullfile(dataDir, 'Wasserstein RDMs');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

addpath('../utils')


%% Load in the neural data

% brain patterns for each subject
data = load(fullfile(dataDir, 'FormattedData_allSubs.mat'));
subs = fieldnames(data.BrainData);
nSubs = length(subs);
fprintf('...loaded brain patterns for all subs (N=%d)!\n', nSubs)
data
data.BrainData

% initial voxel size (isotropic):
voxSize = 3.0; % mm


%% set up the condition comparisons

compPairs = getUniquePairs(data.nConds);
% ...output: pairs x 2 matrix (each row = cond 1 & cond 2 for the pair)
fprintf('...set up pairs to compare %d conditions!\n', data.nConds)



%% loop through the subs

rdmCube = nan(data.nConds, data.nConds, nSubs);
parfor_progress(nSubs);
parfor s = 1:nSubs
   
   %% pull out this sub's data:
   
   currSub = subs{s};
   currData = data.BrainData.(currSub);
   
   rdm = nan(data.nConds, data.nConds);
   
   %% Make meta-voxels
   
   % WD is really computationally-intensive to calculate over large samples,
   % so down-sample the original (3mm x 3mm x 3mm) voxels into "meta-voxels,"
   % where each meta-voxel is 2x2x2 original voxels. Do this oncer per sub (bc
   % different set of voxels for every sub)
   
   [mv_to_v_mat, mv_distmat] = makeMetaVoxels(voxSize, currData.voxelCoords);
   % disp('made meta-voxels.')
   % output:
   % (1) sparse matrix (original voxels x meta-voxels) marking members of each
   % meta-voxel
   % (2) distance matrix for the meta-voxels
   
   % specify the constraints for the model (once per sub):
   dbstop if error
   [A, cTranspose] = makeConstraintMatrix(mv_distmat);
   % disp('made A and c-transpose')
   
   % % eventually: make a sparse version of this
   % [A, cTranspose] = makeSparseConstraintMatrix(mv_distmat, sn);
   % disp('made sparse A and c-transpose')
   
   %% loop through the comparisons

   for c = 1:length(compPairs)
       
       % pull out the activation patterns for the 2 conditions being
       % compared:
       cond1 = compPairs(c, 1);
       cond2 = compPairs(c, 2);
       bp1 = currData.betas(:, cond1);
       bp2 = currData.betas(:, cond2);
       assert(size(bp1, 1) == size(currData.voxelCoords, 1), 'try transposing data matrices.')
       
       %% format the data
        
       disp('formatting the data:')
       
       % normalize:
       [bp1Norm, bp2Norm] = normalize_wasserstein(bp1, bp2);
       % output: 2 vectors (normalized activations for cond 1 and cond 2)
       disp('...normalized')
       
       % downsample:
       dbstop if error
       [bp1MV, bp2MV] = metavoxelize(bp1Norm, bp2Norm, mv_to_v_mat);
       % output: 2 vectors (downsampled activations for cond 1 and cond 2)
       disp('...down-sampled.')
       
       %% wasserstein distance
       wd = computeWasserstein(bp1MV, bp2MV, A, cTranspose);
       % store it twice (for easier checking later on)
       rdm(cond1, cond2) = wd;
       rdm(cond2, cond1) = wd;
       % output: avg. cost (in mm of brain) to move one unit of normalized
       % activation from group 1 to group 2 for this condition.
   end
   
   % store the filled RDM for this sub:
   rdmCube(:, :, s) = rdm;

   parfor_progress; % Count this iteration as progress
   
end

% save the rdm cube:
ConditionPairs = compPairs;
save(fullfile(saveDir, 'wassersteinRDMs-allSubs.mat'), 'rdmCube', 'ConditionPairs');
disp('saved the wasserstein RDMs!')
