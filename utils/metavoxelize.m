function [mvActivations1, mvActivations2] = metavoxelize(activations1, activations2, mv_to_v_mat)
% [mvActivations1, mvActivations2] = metavoxelize(activations1, activations2, mv_to_v_dict, vox_coords)
% Leyla Tarhan & Evan Fields
% 2/2019
% MATLAB R2017b

% downsample activations data from original voxel granularity to meta-voxel 
% granularity by summing across voxels in each meta-voxel.

% inputs:
    % - activations1: betas for group split 1 for 1 condition (# voxels x # 1)
    % - activations2: betas for group split 2 for 1 condition(# voxels x # 1)
    % - mv_to_v_mat: sparse array storing which original voxels are in
    % which meta-voxels (# meta-voxels x # original voxels). 1: original
    % voxel is in the corresponding meta-voxel.
    % - vox_coords: # original voxels x 3 matrix where row i is the x,y,z 
    % coordinates of voxel i for the voxels you're analyzing. In same order
    % as activations1 and activations2
   
% outputs:
    % - mvActivations1: downsampled activations1 (# metavoxels x #
    % 1)
    % - mvActivations2: downsampled activations2
    
 
%% handle inputs

% check the dimensions:
nOrigVoxels = size(mv_to_v_mat, 2);
assert(all(size(activations1) == [nOrigVoxels, 1]), 'WARNING: check orientation of activations matrices.')
assert(all(size(activations2) == [nOrigVoxels, 1]), 'WARNING: check orientation of activations matrices.')

n_mv = size(mv_to_v_mat, 1);
assert(size(mv_to_v_mat, 2) == nOrigVoxels, 'WARNING: unexpected dimensions for matrix of metavoxel memberships.') 

%% downsample betas
% tic
% loop through the groups
for g = 1:2
    fprintf('downsampling activations for group %d...\n', g)
    % handle input -- assumes activations = nVoxels x 1 condition
    activations = [];
    if g == 1
        activations = activations1;
    elseif g == 2
        activations = activations2;
    end
    
    % multiply membership matrix by activations matrix -- because we're
    % aggregating by summing over activations in every original voxel in
    % the meta-voxel, matrix multiplication is a fast way to down-sample
    % (each row of of the membership matrix = 1's for original voxels that
    % are members. Multiply that by the activations = sum(activations that
    % matter*1). This wouldn't be possible if we were aggregating via some
    % other operation, like max or average. 
    mv_acts = mv_to_v_mat*activations;
    
    % handle output
    if g == 1
        mvActivations1 = mv_acts;
    elseif g == 2
        mvActivations2 = mv_acts;
    end
    
end

% fprintf('meta-voxelizing the data took %d minutes.\n', toc/60)

end