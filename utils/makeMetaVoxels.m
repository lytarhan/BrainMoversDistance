function [mv_to_v_mat, mv_distmat] = makeMetaVoxels(VOXEL_SIZE_MM, vox_coords)
% Leyla Tarhan & Evan Fields
% 2/2019
% MATLAB R2017b

% down-sample original voxels into meta-voxels, where each meta-voxel =
% 2x2x2 original voxels. 

% inputs:
    % - VOXEL_SIZE_MM: width of an original voxel, in mm (assuming
    % isotropic)
    % - vox_coords: # voxels x 3 matrix where row i is the x,y,z 
    % coordinates of voxel i for the voxels you're analyzing
    
% - outputs:
    % - mv_t_v_mat: sparse matrix (# metavoxels x # original voxels). 
    % 1: original voxel is a member of corresponding meta-voxel. 
    % - mv_distmat: square distance matrix (in mm) between all pairs of
    % meta-voxels

%% setup

% brain dimensions:
vox_extent = [58 40 46]; % hard-coding: in original voxels
metavox_extent = ceil(vox_extent./2); % in meta-voxels: divide by # of voxels per meta-voxel in each dimension


%% associate original voxels with meta-voxels

% we'll always identify each meta-voxel by its lower-indexed corner
% each metavoxel at [1,1,1] goes from [1,1,1] to [2,2,2]
% since we start indexing at [1,1,1], the lower-indexed corner is also the odd-indexed corner

% only doing this for voxels in our "region of interest" (e.g., reliable
% voxels)

% loop through the voxels and decrease any even coordinates by 1
tic

% set up a dictionary to store the metavoxel coordinates for each original
% voxel:
mv_to_v_dict = containers.Map();


for v = 1:size(vox_coords, 1) % loop through original voxels
    coords = vox_coords(v, :); % pull out this voxel's originak i, j, k coordinates
    metavox_coords = [];
    for i = 1:length(coords) % loop through i, j, k
        if mod(coords(i), 2) == 0 % is it even?
            metavox_coords(i) = coords(i) - 1; % subtract 1 and re-store the coordinate 
        else % it's odd -- keep as-is
            metavox_coords(i) = coords(i);
        end
    end
    
    % store mv and v coordinates as key-value pair
    currKey = num2str(metavox_coords);
    if isKey(mv_to_v_dict, currKey) % already added the key
        % add to the pre-existing values
        mv_to_v_dict(currKey) = [mv_to_v_dict(currKey), {num2str(coords)}];
    else % new key
        mv_to_v_dict(currKey) = {num2str(coords)}; % link to original voxel's coordinates
    end
    
end

% keys = meta-voxel coordinates
% values = original voxel coordinates


%% calculate meta-voxel distance matrix

% get unique metavoxels
n_metavoxels = double(mv_to_v_dict.Count());
mv_coords_unique = keys(mv_to_v_dict);

% set up the distance matrix
mv_distmat = zeros(n_metavoxels, n_metavoxels);

for i = 1:n_metavoxels
    if mod(i, 10) == 0
        fprintf('calculating all distances for meta-voxel %d of %d...\n', i, n_metavoxels);
    end
    for j = (i+1):n_metavoxels
        % get coordinates for voxel i and j:
        mv_coords_i = str2num(mv_coords_unique{i});
        mv_coords_j = str2num(mv_coords_unique{j});
        
        % calculate distance:
        coordsDiff = mv_coords_i - mv_coords_j; 
        dIdx = sqrt(sum(coordsDiff.^2)); % diff in meta-voxel coords
        dMM = dIdx*VOXEL_SIZE_MM; % convert to distance in millimeters
        % multiply by original voxel size bc meta-voxel coordinates are in 
        % same space as original voxels
        
        % store distance:
        mv_distmat(i, j) = dMM;
        mv_distmat(j, i) = dMM;
    end
    
end

% check it out:
% imagesc(mv_distmat); axis square tight; colorbar(); xlabel('metavoxels'); ylabel('metavoxels')

%% store meta-voxel memberships

% take the dictionary and make it into a sparse matrix storing the original
% voxels that are members of the meta-voxels:
disp('...storing meta-voxel memberships in a sparse matrix...')
n_origvoxels = size(vox_coords, 1);
mv_to_v_mat = sparse(n_metavoxels, n_origvoxels); % original voxels x meta-voxels

% loop through the meta-voxels
for mv = 1:n_metavoxels
   % get current mv coordinates:
   currKey = mv_coords_unique{mv};
   mv_coords = str2num(currKey); % pull back into number format
   
   % get the original voxels in this meta-voxel:
   currValues = mv_to_v_dict(currKey);
   % figure out the index corresponding to those voxels:
   for v = 1:length(currValues)
      % match by coordinates
      ov_coords = str2num(currValues{v});
      ovIdx = find(ismember(vox_coords, ov_coords, 'rows'));
      
      % update the membership matrix:
      mv_to_v_mat(mv, ovIdx) = 1;
       
   end
    
end


fprintf('making %d voxels into %d meta-voxelsde meta-voxels took %0.2f s\n', size(vox_coords, 1), n_metavoxels, toc)

end
