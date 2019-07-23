function wd = computeWasserstein(activations1, activations2, A, cTranspose)
% Leyla Tarhan & Evan Fields
% 2/2019
% MATLAB R2017b

% calculate wasserstein (aka earth-mover's) distance between two activation
% patterns.

% to use this code, you'll need to have Gurobi installed  with an academic 
% license and add its directory to MATLAB's path. 
% For more info on this, see http://www.gurobi.com/documentation/

% inputs:
    % - activations1 & activations2: vectors of betas in all the voxels 
    % you're analyzing for a single condition -- 1 per group. These
    % patterns should have the same "mass" (sum) and only non-negative
    % values.
    % - A: constraint matrix, made using makeConstraintMatrix() or
    % makeSparseConstraintMatrix()
    % - cTranspose: c-transpose vector, made using makeConstraintMatrix()
    % or makeSparseConstraintMatrix()
    
% output:
    % - wd: wasserstein distance between the 2 groups' activation patterns.
    % Technically the optimal cost of transforming one 3-D distribution of
    % data into the other. Units = avg. mm of brain (or whatever your distance
    % metric is) needed to transform one unit of activation (which might be
    % normalized) from group 1 to group 2.
    
%% Handle the inputs

% make sure activations1 and activations2 have roughly equal mass
tolerance = 1e-7; % reasonable margin of error
diffData = sum(activations1) - sum(activations2);
assert(abs(diffData)<=tolerance, 'brain patterns have unequal mass')

% figure out which data vector is "lighter" (calculate the cost to transform this one to
% the other)
if sum(activations1) < sum(activations2) % activations1 is lighter
    x = activations2; % always make x the heavier vector
    y = activations1;
else
    x = activations1;
    y = activations2;
end

% what if sums are actually equal? Then it doesn't matter which is x and
% which is y

n = length(x);
assert(n == length(y), 'x and y are different lengths')

%% set up the model

% % debugging mode -- cut the data down to a smaller problem:
% disp('using a toy dataset.')
% nVox_debug = 1500;
% x = x(1:nVox_debug);
% y = y(1:nVox_debug);
% % normalize so these add to 1 (since just lifted them from the normalized
% % data)
% x = x./sum(x);
% y = y./sum(y);
% n = length(x);

% specify the model
% A: matrix of coefficients on the problem's constraints (# metavoxels x #
% ordered pairs of metavoxels)
% x: vector of decision variables for each ordered pair of meta-voxels (aka
% the amount of flow from meta-voxel i to meta-voxel j)
% - A*x = b, where b = a vector with (group y activation - group x
% activation) for each metavoxel.
% - rows of A sum to 0
% - values in A = 0, -1, or 1. Each element modifies the activation that 
% flows from meta-voxel i to meta-voxel j. So, if activation should
% flow out of mv i into mv j, a for [i, j] will = 1 and a for [j, i] will =
% -1. 
model = [];
model.A = A; % constraint matrix (has to be a sparse matrix)
assert(size(A, 1) == n, 'constraint matrix doesn''t match data size.')
model.obj = cTranspose; % objective function (what we're trying to minimize)
% actually specifying a vector of values (c-transpose), aka the distances 
% between each ordered pair of meta-voxels. 
% The objective function is c-transpose * x
assert(sum(isnan(model.obj)) == 0, 'c-transpose vector skipped some elements.')

% b vector (difference in activations between groups)
model.rhs = x-y; 
% modeling flow from x to y.
% if x & y really have the same mass, this really only matters for
% readability:
    % - % In A, each row = net *outflow* from each voxel.
    % - here, net outflow = x activation - y activation

% other attributes:
model.sense = '='; % constraint sense vector (whether each constraint in A
% is an equality or inequality constraint
model.modelsense = '1'; % sense of the objective function: minimization
model.vtype = 'C'; % continuous variables
model.LB = 0; % lower bound on amount of flow moved between voxels
% http://www.gurobi.com/documentation/8.1/refman/parameters.html#sec:Parameters
% for more

%% solve the model

params.outputflag = 1;
% tic
disp('solving the model...')
result = gurobi(model, params);
% fprintf('Solved the model with %d meta-voxels in %0.2f s!\n', n, toc)

% extract the output
wd = result.objval;


end

