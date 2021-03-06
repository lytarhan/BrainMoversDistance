% Leyla Tarhan
% https://github.com/lytarhan
% 10/2019
% MATLAB R2017b

% Step 1 in use case #1 (replicability analysis): make splits of the
% example data. 

% generate 100 random splits among the subjects in the data. On each iteration,
% randomly select 6 subs to be in group 1 and 7 subs to be in group 2. Then
% average across all the subjects in each group and save the result as a 
% .mat file. 

% To use this script: download the dataset "1-Replicability" from this project's OSF
% repository (https://osf.io/m9ac3/) to a directory called "Data", which 
% should be saved to the same directory as this script.

%% clean up

clear all
close all
clc

%% set up the file structure

% single-sub betas:
dataDir = 'Data';

% where to save the split data:
saveDir = fullfile(dataDir, 'split averages');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% which subs?
subs = [1:13];
nSubs = length(subs);

% how many splits?
nSplits = 100;

% group sizes
group1N = 6;
group2N = 7;
assert((group1N + group2N) == nSubs, 'subs vector size doesn''t match specified group sizes.')

% how many voxels?
nVoxels = 8125;

% how many conditions in the experiment (videos)?
nConds = 120;

%% load the betas

betasCube = zeros(nVoxels, nConds, nSubs);
for s = 1:nSubs
    subNum = subs(s);
    subName = sprintf('Sub%d', subNum);
    
    % load their betas:
    data = load(fullfile(dataDir, ['SSBetas_' subName '.mat']));
    
    % concatenate set 1 and set 2 data:
    betasAllConds = [data.B.set1.betas, data.B.set2.betas];
    assert(all(size(betasAllConds) == [nVoxels, nConds]), 'betasAllConds dimensions are off.')
    
    % add to the cube:
    betasCube(:, :, s) = betasAllConds;
    
end

disp('loaded betas for all subs.')

%% make random splits

disp('computing group splits...')

for i = 1:nSplits
    fprintf('computing split %d of %d...\n', i, nSplits)
    group1Avg = [];
    group2Avg = [];
    
    % randomly select subs:
    randOrder = randperm(nSubs);
    
    % group 1 = first group1N
    group1Subs = randOrder(1:group1N);
    assert(length(group1Subs) == group1N, 'group1 is the wrong size.')
    group1Cube = betasCube(:, :, group1Subs);
    % group 1 average
    group1Avg = mean(group1Cube, 3);
    
    % group 2 = remaining group2N
    group2Subs = randOrder(group1N+1:end);
    assert(length(group2Subs) == group2N, 'group2 is the wrong size.')
    assert(isempty(intersect(group1Subs, group2Subs)), 'group 1 and group 2 have overlapping subs!')
    group2Cube = betasCube(:, :, group2Subs);
    % group 2 average
    group2Avg = mean(group2Cube, 3);
    
    % save it:
    fn = sprintf('groupAverages_split%d.mat', i);
    save(fullfile(saveDir, fn), 'group1Avg', 'group2Avg')
    
end

disp('...computed all splits!')