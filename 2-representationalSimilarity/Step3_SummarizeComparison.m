% Leyla Tarhan
% https://github.com/lytarhan
% 1/2020
% MATLAB R2017b


% Step 3 in use case #2 (representational similarity analysis): 
% compare brain-behavior correlations when neural dissimilarity is computed
% using correlation distance vs. wasserstein distance. Compare them with 
% reference to the results from Step 2b, when the matrices were scrambled. 
% Put another way, which distance metric produces a neural dissimilarity 
% matrix that more closely resembles the behavioral dissimilarity matrix, 
% compared to a scrambled baseline?


%% clean up

clear all
close all
clc

%% file structure

% results from RSA:
dataDir = 'Results';

% developing over-ride:
devFlag = 1;
if devFlag
    dataDir = 'C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Project - BrainMoversDistance\Outputs\OSF - DataForGitHub\2-RepresentationalSimilarity\Results';
end


addpath('../utils')



%% load in RSA results
% from Steps 2a and 2b

% intact data:
intactData = load(fullfile(dataDir, 'Brain-VisSearchRSA.mat'));

% scrambled data:
scramData = load(fullfile(dataDir, 'Brain-VisSearchRSA_Scrambled.mat'));

%% Summarize results

% (1) check out all 4 bars at once, to get a sense of how the intact vs. scrambled data compare:
intactScramData = [mean(intactData.rsaCorr), mean(mean(scramData.rsaScramCorr)); mean(intactData.rsaWD), mean(mean(scramData.rsaScramWD))];
b1 = bar(intactScramData);
legend({'intact', 'scrambled'});
ylabel('brain-behavior correlation');
xticklabels({'correlation distance', 'wasserstein distance'})

% (2) plot intact - scrambled for each sub:
% correlation-distance
rDiffVsScram = intactData.rsaCorr - mean(scramData.rsaScramCorr, 2); % average over scrambling iterations
% wasserstein distance
wdDiffVsScram = intactData.rsaWD - mean(scramData.rsaScramWD, 2); 

% plot it:
figure('Color', [1 1 1])
b2 = boxplot([rDiffVsScram, wdDiffVsScram]); % correlation, the wasserstein
ylabel('brain-behavior correlation (intact - scrambled)')
xlabel('neural dissimilarity measure')
xticklabels({'correlation distance', 'wasserstein distance'})

% save it:
saveFigureHelper(1, dataDir, 'Brain-VisSearch_corrVSwasserstein_vsScram.png')

%% stats
% Q: compared to scrambled-model baseline, is the brain-behavior 
% correlation higher when you measure neural dissimilarity using 
% correlation distance or wasserstein distance?
clc
[h, p, ci, stats] = ttest(rDiffVsScram, wdDiffVsScram);
if h > 0
    disp('there''s a difference between correlation and wasserstein distance, vs. scrambled baseline!');
    if stats.tstat > 0 % ttest(x, y): t-stat is based on x - y
        disp('correlation is higher');
    else
        disp('wasserstein is higher');
    end
    fprintf('t(%d) = %0.2f, p = %0.4f\n', stats.df, stats.tstat, p)
else
   disp('there''s no difference between correlation and wasserstein distance.') 
end


