function [A, cTranspose] = makeConstraintMatrix(distMat, speedyFlag)
% Leyla Tarhan & Evan Fields
% 2/2019
% MATLAB R2017b

% compute the constraint matrix and c-transpose array to feed into a solver
% later on while solving wasserstein distances. Only have to make these
% once per analysis (after making the meta-voxels).

% NB: takes awhile (on the order of 20 mins)!

% inputs:
    % - distMat: square distance matrix for all the voxels you're
    % analyzing, in actual distance (e.g., millimeters).
    % - speedyFlag: use speedier method to make the A matrix? default = 1
    
% outputs:
    % - % A: matrix of coefficients on the problem's constraints (# metavoxels x #
    % ordered pairs of metavoxels)
    % - cTranspose: a vector of values (c-transpose), aka the distances 
    % between each ordered pair of meta-voxels. Objective function we're
    % eventually trying to minimize is c-transpose * flow in/out of all
    % voxels
    
%% setup

if ~exist('speedyFlag', 'var')
    speedyFlag = 1;
end

% make sure distMat is square
assert(size(distMat, 1) == size(distMat, 2), 'distMat should be square.')

% # meta-voxels
n = size(distMat, 1);

% % debugging version:
% disp('toy dataset.')
% n = 1500;

% constraint matrix (has to be a sparse matrix):
if speedyFlag
    % set up to fill in a sparse matrix -- much faster to store locations and
    % values of non-zero elements than to store/update a sparse matrix in the
    % loop below
    vals = []; % array to store constraint values
    rowIdx = []; % array to store row indices for constraint values
    colIdx = []; % array to store column indices for constraint values
else % normal sparse matrix
    A = sparse(n, n^2);
end
    

cTranspose = NaN(1, n^2); % objective function (what we're trying to minimize)
% actually specifying a vector of values (c-transpose), aka the distances 
% between each ordered pair of meta-voxels. 
% The objective function is c-transpose * x

%% Build up A and c-transpose

% take the lower triangle of our square distance matrix, *including* the 
% main diagonal, and unfold it into a vector. Unfold it by hand rather than
% using squareform() to ensure that the order is the same across c and x.

disp('Specifying constraint matrix and c-transpose...')
tic
for i = 1:n
   for j = 1:n
      % fill in c-transpose vector with the distance between these voxels
      currIdx = ijToIdx(i, j, n);
      cTranspose(currIdx) = distMat(i, j);
      
      % fill in A matrix
      if i == j % flow is leaving and arriving at the same voxel
          continue % make sure not to overwrite values (a stays = 0, as initialized)
      else
          % flow from i --> j
          idx_iToj = ijToIdx(i, j, n);
          if speedyFlag
              rowIdx = [rowIdx, i];
              colIdx = [colIdx, idx_iToj];
              vals = [vals, 1]; % a = 1
          else
              A(i, idx_iToj) = 1; % a = 1
          end
          
          % flow from j --> i
          idx_jToi = ijToIdx(j, i, n);
          if speedyFlag
              rowIdx = [rowIdx, i];
              colIdx = [colIdx, idx_jToi];
              vals = [vals, -1]; % a = -1
          else
              A(i, idx_jToi) = -1; % a = -1
          end
          % each row's dot product = net *outflow* from the corresponding
          % voxel
      end     
      
   end
end
if speedyFlag
    % make the sparse A matrix:
    A = sparse(rowIdx, colIdx, vals, n, n^2);
end
fprintf('constraint matrix for %d voxels took %0.2f s\n', n, toc)



end

% -------------------------------------------------------------------------

%% 
function idx = ijToIdx(i, j, nVox)
% given (i, j) return the index in x that corresponds to the flow from i to
% j
idx = nVox*(i-1) + j; % make vector by stacking rows of the square matrix
end

function [i, j] = idxToij(idx, nVox)
% given idx return (i, j) which correspond to the meta-voxels whose flow is
% specified by x(idx).
j = mod(idx, nVox);
i = (idx - j)/nVox + 1;
end