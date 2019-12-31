function [A, cTranspose] = makeSparseConstraintMatrix(distMat, sparseNet)
% Leyla Tarhan & Evan Fields
% 2/2019
% MATLAB R2017b

% compute the sparse constraint matrix and c-transpose array to feed into a solver
% later on while solving wasserstein distances. Only have to make these
% once per analysis (after making the meta-voxels).

% inputs:
    % - distMat: square distance matrix for all the voxels you're
    % analyzing, in actual distance (e.g., millimeters).
    % - sparseNet: matrix (2 x # edges in the sparse network). Each row =
    % an edge between metavoxel in col. 1 and that in col. 2.

% outputs:
    % - % A: matrix of coefficients on the problem's constraints (# metavoxels x #
    % ordered pairs of metavoxels). Because this is a sparse version, many 
    % columns of A will = 0 (bc these edges won't be used in the wasserstein
    % computation).
    % - cTranspose: a vector of values (c-transpose), aka the distances 
    % between each ordered pair of meta-voxels. Objective function we're
    % eventually trying to minimize is c-transpose * flow in/out of all
    % voxels. Because this is a sparse version, many columns of c_Transpose
    % will = 0 (bc these edges won't be used in the wasserstein
    % computation).
    
%% setup

% convert sparse network to array:
sparseNet = table2array(sparseNet);

% make sure distMat is square
assert(size(distMat, 1) == size(distMat, 2), 'distMat should be square.')

% # meta-voxels
n = size(distMat, 1);
% # edges between meta-voxels in the sparse network
nEdges = size(sparseNet, 1);

% debugging version:
% disp('toy dataset.')
% n = 100;

% constraint matrix (has to be a sparse matrix):

% set up to fill in a sparse matrix -- much faster to store locations and
% values of non-zero elements than to store/update a sparse matrix in the
% loop below
vals = NaN(nEdges*4, 1); % array to store constraint values
rowIdx = NaN(nEdges*4, 1); % array to store row indices for constraint values
colIdx = NaN(nEdges*4, 1); % array to store column indices for constraint values
    
cTranspose = zeros(1, n^2); % objective function (what we're trying to minimize)
% actually specifying a vector of values (c-transpose), aka the distances 
% between each ordered pair of meta-voxels. 
% The objective function is c-transpose * x

%% Build up A and c-transpose

% take the lower triangle of our square distance matrix, *including* the 
% main diagonal, and unfold it into a vector. Unfold it by hand rather than
% using squareform() to ensure that the order is the same across c and x.

disp('Specifying constraint matrix and c-transpose...')
tic
% loop through the sparse network
constraintIdx = 1;
for e = 1:nEdges
    i = sparseNet(e, 1);
    j = sparseNet(e, 2);
    
    % fill in c-transpose vector with the distance between these voxels
    idx_iToj = ijToIdx(i, j, n); % index for flow from i --> j
    cTranspose(idx_iToj) = distMat(i, j);
    idx_jToi = ijToIdx(j, i, n); % index for flow from j --> i
    cTranspose(idx_jToi) = distMat(j, i);
    
    % fill in A matrix
    if i == j % flow is leaving and arriving at the same voxel
       continue
    else
        % out-flow for voxel i
        rowIdx(constraintIdx) = i;
        colIdx(constraintIdx) = idx_iToj;
        vals(constraintIdx) = 1; % a = 1
        constraintIdx = constraintIdx + 1;
        
        % in-flow for voxel i
        rowIdx(constraintIdx) = i;
        colIdx(constraintIdx) = idx_jToi;
        vals(constraintIdx) = -1; % a = -1
        constraintIdx = constraintIdx + 1;
        
        % out-flow for voxel j
        rowIdx(constraintIdx) = j;
        colIdx(constraintIdx) = idx_jToi;
        vals(constraintIdx) = 1; % a = 1
        constraintIdx = constraintIdx + 1;
        
        % in-flow for voxel j
        rowIdx(constraintIdx) = j;
        colIdx(constraintIdx) = idx_iToj;
        vals(constraintIdx) = -1; % a = -1
        constraintIdx = constraintIdx + 1;

        % each row's dot product = net *outflow* from the corresponding
        % voxel
    end
end

% make the sparse A matrix:
A = sparse(rowIdx, colIdx, vals, n, n^2);

fprintf('sparse constraint matrix for %d voxels took %0.2f s\n', n, toc)



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