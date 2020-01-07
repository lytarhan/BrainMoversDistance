% Leyla Tarhan
% https://github.com/lytarhan
% 1/2020
% MATLAB R2017b


% Step 2b in use case #2 (representational similarity analysis): relate 
% neural and behavioral dissimilarities (separately for each subject, in a 
% pre-defined brain region of interest). Do this with neural 
% dissimilarities calculated using correlation distance and Wasserstein 
% Distance. This time, scrambled the neural dissimilarities and relate them
% to the intact behavioral dissimilarities, to set a baseline for the RSA
% results. 


%% clean up

clear all
close all
clc

%% file structure

topDir = pwd;

bDataDir = fullfile(topDir, 'Data-Behavior'); % behavioral data
wdDataDir = fullfile(topDir, 'Data-fMRI', 'Wasserstein RDMs'); % neural RDMs, made using Wasserstein distance
rDataDir = fullfile(topDir, 'Data-fMRI'); % neural RDMs, made using correlation distance

saveDir = fullfile(topDir, 'Results');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

addpath('../utils')

%% setup

% how many conditions in the data?
nConds = 72;
nPairs = factorial(nConds) / (factorial(2) * factorial(nConds - 2));

% how many iterations for scrambling?
iters = 100;


%% format the data
% output:
% ...group-avg. visual search RDM (higher values = more dissimilar)
% ...correlation RDM for each fMRI sub (1-corr --> higher values = more
% dissimilar)
% ...WD RDM for each fMRI sub (higher values = more dissimlar)


% behavioral data (visual search times):
vsData = load(fullfile(bDataDir, 'SearchData'));
vsRDM = vsData.SearchData.condPairCoeffs; % lower triangle, minus the main diagonal

% neural RDM: correlation distance (1 per sub)
bpData = load(fullfile(rDataDir, 'FormattedData_allSubs.mat'));
subs = fieldnames(bpData.BrainData);
% set up a matrix to hold the RDMs
rDistMat = nan(nPairs, length(subs));
for s = 1:length(subs)
   rDistMat(:, s) = bpData.BrainData.(subs{s}).corrRDM; 
end

% neural RDM: wasserstein distance (1 per sub)
wdData = load(fullfile(wdDataDir, 'wassersteinRDMs-allSubs.mat'));
wdDistMat = nan(nPairs, length(subs));
for s = 1:length(subs)
    wdDistMat(:, s) = getLowerTri(wdData.rdmCube(:, :, s));
end

%% RSA with scrambling: loop through the subs

% set up to store the scrambled-brain-behavior correlations:
rsaScramCorr = nan(length(subs), iters);
rsaScramWD = nan(length(subs), iters);
for s = 1:length(subs)
    % loop through scrambling iterations
    for i = 1:iters
        randIdx = randperm(nPairs); % random scramble -- apply to both kinds of neural RDM
        % correlation distance:
        rsaScramCorr(s, i) = corr(rDistMat(randIdx, s), vsRDM);
        
        % wasserstein distance:
        rsaScramWD(s, i) = corr(wdDistMat(randIdx, s), vsRDM);
    end
end
disp('Got scrambled brain-behavior correlations for all subs!')
rsaScramCorr = rsaScramCorr.*(-1);
rsaScramWD = rsaScramWD.*(-1);
disp('Multiplied correlations by -1 so more interpretable.')


% save it
SearchDataNote = vsData.SearchData.meta;
CorrelationsNote = 'multiplied correlations by -1 bc high values = more similar in search data, but less similar in brains';
save(fullfile(saveDir, 'Brain-VisSearchRSA_Scrambled.mat'), 'rsaScramCorr', 'rsaScramWD', 'subs', 'SearchDataNote', 'CorrelationsNote');
disp('saved scrambled RSA results!')

%% Quick plots
% sanity-check

% average over iterations:
rsaScramCorrAvg = mean(rsaScramCorr, 2);
rsaScramWDAvg = mean(rsaScramWD, 2);

figure('Color', [1 1 1])
b = boxplot([rsaScramCorrAvg, rsaScramWDAvg]);
ylabel('scrambled brain-behavior correlation')
xlabel('neural dissimilarity measure')
xticklabels({'correlation distance', 'wasserstein distance'})

% save it:
saveFigureHelper(1, saveDir, 'Brain-VisSearch_corrVSwd_Scrambled.png')
