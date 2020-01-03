% Leyla Tarhan
% 5/2019
% MATLAB R2017b

% relate scrambled neural and intact behavioral dissimilarities (set a
% baseline for comparing RSA using wasserstein and correlation distance).


% To Do:


%% clean up

clear all
close all
clc

%% setup

% how many conds?
nConds = 72;
nPairs = factorial(nConds) / (factorial(2) * factorial(nConds - 2));

% how many iterations for scrambling?
iters = 100;

%% file structure

topDataDir = 'C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Project - BrainMoversDistance\Experiment - ExploringObjectsRSA\Analysis';

% behavioral data
bDataDir = fullfile(topDataDir, 'Data - Behavior');

% wasserstein RDMs
wdDataDir = fullfile(topDataDir, 'Data - fMRI', 'Wasserstein RDMs');

% correlation-distance RDMs
rDataDir = fullfile(topDataDir, 'Data - fMRI', 'FormattedData');

% figure- and results-saving
saveDir = fullfile(topDataDir, 'Results');
if~exist(saveDir, 'dir'); mkdir(saveDir); end

% helpers:
addpath(genpath('C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Code - frequent helpers'));


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

% make sure behavioral and neural RDMs are in the same order
for p = 1:nPairs
   % which items in the visual search data?
   vs1 = vsData.SearchData.labelIm1{p};
   vs2 = vsData.SearchData.labelIm2{p};
   
   % which items in the neural data?
   currConds = wdData.ConditionPairs(p, :);
   neural1 = bpData.CondNames{currConds(1)};
   neural2 = bpData.CondNames{currConds(2)};
   
   % do they match?
   match1 = strcmp(vs1, neural1) || strcmp(vs1, neural2);
   match2 = strcmp(vs2, neural1) || strcmp(vs2, neural2);
   if ~match1 || ~match2
       disp('mismatch!')
       keyboard
   end
end
disp('neural and behavioral RDMs are in the same order!')

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
