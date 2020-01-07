% Leyla Tarhan
% https://github.com/lytarhan
% 1/2020
% MATLAB R2017b


% Step 2a in use case #2 (representational similarity analysis): relate 
% neural and behavioral dissimilarities (separately for each subject, in a 
% pre-defined brain region of interest). Do this with neural 
% dissimilarities calculated using correlation distance and Wasserstein 
% Distance. 


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


%% compare the RDMs

% quick check: how similar are correlation-based and wasserstein-based
% dissimilarity? Correlate these two RDMs for each subject.
rdmCorrs = nan(length(subs), 1);
for s = 1:length(subs)
   rRDM = rDistMat(:, s);
   wdRDM = wdDistMat(:, s);
   rdmCorrs(s) = corr(rRDM, wdRDM);
   
   % make a quick figure:
   figure('Color', [1 1 1], 'Position', [10 60 800 400])
   subplot(1, 3, 1); imagesc(squareform(rRDM)); colorbar(); axis square tight off; title('correlation RDM');
   subplot(1, 3, 3); imagesc(squareform(wdRDM)); colorbar(); axis square tight off; title('wasserstein RDM');
   subplot(1, 3, 2); h = text(.5, .5, sprintf('<-- r = %0.2f -->', rdmCorrs(s))); axis square tight off
   set(h, 'HorizontalAlignment', 'center');
   % save it:
   sn = ['compareNeuralRDMs_' subs{s} '.png'];
   saveFigureHelper(1, saveDir, sn);
   close
end

figure('Color', [1 1 1])
boxplot(rdmCorrs);
title('comparing wasserstein and correlation RDMs')
ylabel('rdm-to-rdm correlation');
disp('average correlation between correlation-based and wasserstein-based neural RDMs:')
fprintf('%0.2f\n', mean(rdmCorrs))


%% RSA: loop through the subs

% set up to store the brain-behavior correlations:
rsaCorr = nan(length(subs), 1);
rsaWD = nan(length(subs), 1);
for s = 1:length(subs)
    % brain-behavior correlation for corr-distance neural RDM:
    rsaCorr(s) = corr(rDistMat(:, s), vsRDM);
    
    % brain-behavior correlation for WD-distance neural RDM:
    rsaWD(s) = corr(wdDistMat(:, s), vsRDM);
end
disp('Got brain-behavior correlations for all subs!')
rsaCorr = rsaCorr.*(-1);
rsaWD = rsaWD.*(-1);
disp('Multiplied correlations by -1 so more interpretable.')

% save it:
SearchDataNote = vsData.SearchData.meta;
CorrelationsNote = 'multiplied correlations by -1 bc high values = more similar in search data, but less similar in brains';
save(fullfile(saveDir, 'Brain-VisSearchRSA.mat'), 'rsaCorr', 'rsaWD', 'subs', 'SearchDataNote', 'CorrelationsNote');
disp('saved results!')

%% Plot it

figure('Color', [1 1 1])
b = boxplot([rsaCorr, rsaWD]);
ylabel('brain-behavior correlation')
xlabel('neural dissimilarity measure')
xticklabels({'correlation distance', 'wasserstein distance'})

% save it:
saveFigureHelper(1, saveDir, 'Brain-VisSearch_corrVSwasserstein.png')

%% Preliminary stats

% Q: is there a difference in the brain-behavior correlation between the
% two neural dissimilarity metrics?
clc
[h, p, ci, stats] = ttest(rsaCorr, rsaWD);
if h > 0
    disp('there''s a difference between correlation and wasserstein distance!');
    if stats.tstat > 0 % ttest(x, y): t-stat is based on x - y
        disp('correlation is higher');
    else
        disp('wasserstein is higher');
    end
    fprintf('t(%d) = %0.2f, p = %0.4f\n', stats.df, stats.tstat, p)
else
   disp('there''s no difference between correlation and wasserstein distance.') 
end

