% Leyla Tarhan
% 5/2019
% MATLAB R2017b

% compare brain-behavior correlation when neural dissimilarity is computed
% using correlation distance vs. wasserstein distance (compared to a
% scrambled model baseline, which type of neural dissimilarity is more
% similar to behavioral dissimilarity?)

% To Do:


%% clean up

clear all
close all
clc

%% file structure

% results from RSA:
dataDir = 'C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Project - BrainMoversDistance\Experiment - ExploringObjectsRSA\Analysis\Results';

% helpers:
addpath(genpath('C:\Users\Leyla\Dropbox (KonkLab)\Research-Tarhan\Code - frequent helpers'));

%% load in RSA results

% intact data:
intactData = load(fullfile(dataDir, 'Brain-VisSearchRSA.mat'));

% scrambled data:
scramData = load(fullfile(dataDir, 'Brain-VisSearchRSA_Scrambled.mat'));

%% Summarize results

% (1) check out all 4 bars at once:
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


