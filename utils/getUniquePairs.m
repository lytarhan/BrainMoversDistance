function pairsMat = getUniquePairs(nConds)
% Leyla Tarhan
% 4/2019
% MATLAB R2017b

% given n conditions, figure out all the unique, non-redundant pairs of 2
% conditions. 

% output: matrix of condition #'s (# pairs x 2): each row = cond 1 & cond 2 for the pair

%%
% start with a cube (conds x conds x order in the pair):
compsCube = zeros(nConds, nConds, 2);
for i = 1:nConds
    for j = 1:nConds
       compsCube(i, j, 1) = i; % condition 1 being compared
       compsCube(i, j, 2) = j; % condition 2 being compared
    end   
end

% % check it out:
% figure('Position', [10 60 800 400])
% subplot(1, 2, 1); imagesc(compsCube(:, :, 1)); axis square; title('position 1'); colorbar()
% subplot(1, 2, 2); imagesc(compsCube(:, :, 2)); axis square; title('position 2'); colorbar()

% get the unique, non-redundant pairs: grab just the lower triangle, minus the main diagonal
pairsMat = [];
for x = 1:2
    pairsMat(:, x) = getLowerTri(compsCube(:, :, x));
end

% % check it out:
% imagesc(pairsMat); colorbar(); ylabel('pair'); title('setting up condition pairs')
% xticks(1:2); xticklabels({'position 1', 'position 2'})

end