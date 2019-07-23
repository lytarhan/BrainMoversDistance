% Leyla Tarhan
% https://github.com/lytarhan
% 8/2019
% MATLAB R2017b

% Step 4 in use case #1 (replicability analysis): check out the results of 
% the analysis!

% plot wasserstein distance for each condition, averaged across 100 splits.
% Then compare to the scrambled results' baseline to determine whether the 
% data were replicable. 

%% clean up

clear all
close all
clc

%% folder structure

% directory with original data:
dataDir = 'Data';

% directories with wasserstein output:
resultsDir = 'Results'; % intact data
resultsScramDir = 'Results-Scrambled'; % scrambled data

% directory to save figures in:
saveDir = 'Figures';
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% helpers:
addpath('../utils');

%% Load in the data

% condition names:
load(fullfile(dataDir, 'SSBetas_Sub3.mat')); % example sub for condition names
condNames = [B.set1.condNames, B.set2.condNames];
condNamesClean = [];
for c = 1:length(condNames)
    condNameSplit = strsplit(condNames{c}, '_');
    condNamesClean{c} = condNameSplit{2};
end
disp('loaded condition names!')


% wasserstein results:
nSplits = 100;
nConds = 120;

wdMat = zeros(nConds, nSplits);
wdScramMat = zeros(nConds, nSplits);

for i = 1:nSplits
   fn = sprintf('wasserstein_split%d.mat', i);
   
   % intact results:
   wd = load(fullfile(resultsDir, fn));
   wdMat(:, i) = wd.wd;
   
   % scrambled results:
   wdScram = load(fullfile(resultsScramDir, fn));
   wdScramMat(:, i) = wdScram.wd;
   
   clear wd wdScram
end
fprintf('loaded intact and scrambled distances for all %d splits!\n', nSplits)

%% intact data: plot mean by condition

% average:
wdAvg = mean(wdMat, 2);

% sort it: 
[s, si] = sort(wdAvg, 'descend');

% plot it:
figure('Color', [1 1 1], 'Position', [10 60 1350 400]);
b = bar(s);
title(sprintf('Wasserstein Distances: %d splits', nSplits));
ylabel('Avg. Distance (mm)')
xticks(1:nConds); xticklabels(condNamesClean(si)); xtickangle(50) 
ax = ancestor(b, 'axes');
ax.XAxis.FontSize = 8;

% add error bars:
hold on
se = std(wdMat')./sqrt(nSplits); % standard error
eb = errorbar(s, se(si)'); % sorted 
hold off

% save the figure:
sn = sprintf('WDByCondition_%dSplits.png', nSplits);
saveFigureHelper(1, saveDir, sn);
close

% print out the grand mean:
clc
fprintf('Overall average wasserstein distance: %0.2f mm\n', mean(wdAvg))
disp('units: average mm to move 1 unit of normalized activation from group 1 to group 2.')


%% Interpretation: compare intact and scrambled results

wdAvg = mean(wdMat, 2);
wdScramAvg = mean(wdScramMat, 2);

figure('Color', [1 1 1], 'Position', [10 60 900 400])
% means
subplot(1, 2, 1)
b = bar([mean(wdAvg), mean(wdScramAvg)]);
title('Wasserstein')
ylabel('Wasserstein Distance (mm)');
xticklabels({'intact', 'scrambled'});
ylim([0, 11]);
b.FaceColor = 'flat';
b.CData(1,:) = [0 100/255 0]; % green for WD
b.CData(2,:) = [180 180 180]./255; % grey for scrambled

% error bars
wdSE = std(wdAvg)./sqrt(nConds); % standard error
wdScramSE = std(wdScramAvg)./sqrt(nConds);
hold on
eb = errorbar([mean(wdAvg), mean(wdScramAvg)], [wdSE, wdScramSE]); 
eb(1).LineStyle = 'none';
hold off



% stats (compare distances from intact vs. scrambled data)
wdScramAvg = mean(wdScramMat, 2);
wdAvg = mean(wdMat, 2);

% paired t-test (pair the conditions)
[H, P, CI, stats] = ttest(wdScramAvg, wdAvg); 

% print it out:
fprintf('intact vs. scrambled results:\n')
fprintf('===================================================\n')
fprintf('t(%d) = %0.2f, p = %0.3f \n', stats.df, stats.tstat, P)
if P < 0.001
    disp('* test is significant at p<0.001')
end


