function [group1Norm, group2Norm] = normalize_wasserstein(group1, group2)
% Leyla Tarhan & Evan Fields
% 2/2019
% MATLAB R2017b

% normalize 2 data vectors in order to compute wasserstein distance between
% them.

% inputs:
    % - group1 & group2: 2 data vectors (activations in each voxel to a
    % single condition)
    
% outputs: 
    % - group1Norm & group2Norm: 2 normalized data vectors (vectors have
    % equal "mass" and only values >= 0)
    
%% 
% tic
% normalize vec 1
group1Norm = group1 - min(group1); % shift everything to be > 0
group1Norm = group1Norm./(sum(group1Norm)); % normalize to range [0, 1] so equal masses

% normalize vec 2
group2Norm = group2 - min(group2); % shift everything to be > 0
group2Norm = group2Norm./(sum(group2Norm)); % normalize to range [0, 1] so equal masses

% check that it worked:
assert(all(group1Norm >=0) && all(group2Norm >= 0), 'detected some negative values.')
% assert(sum(group1Norm) == sum(group2Norm), 'detected unequal masses.')

% fprintf('data normalization took %d seconds.\n', toc)
end