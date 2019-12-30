function [updatedNet, updatedWeightedNet] = trimEdges(sparseNet, voxDists, allowed_inflation)
% Leyla Tarhan & Evan Fields
% trim unnecessary edges from the sparse network (generally follows the
% edge-adding step). 

% inputs:
    % ...sparseNet: nVoxels x nVoxels. If 2 voxels are directly connected
    % by an edge, the corresponding cell = 1.
    % ...voxDists: nVoxels x nVoxels. Each cell = actual distance (in the
    % brain, mm) between voxels i & j.
    % ...allowed_inflation: parameter to decide when to trim edges. Should
    % be > 1 but not by much.

% output:
    % ...updatedNet: updated version of sparseNet, with fewer edges.

%% setup

fprintf('\ntrimming edges (be patient, this will take awhile)...\n')

% how many voxels?
nVox = size(sparseNet, 1);

% setup to store the edges being trimmed (should be faster than updating
% the sparse net directly, even using some tricks that still require some
% really long vectors):
trimCells = [];

% get the locations of the existing net's edges:
[prevRow, prevCol] = find(sparseNet); 

% flag for shorter version to help in developing code:
devFlag = 1;
if devFlag
    disp(' ')
    disp('RUNNING IN DEV MODE.')
    disp(' ')
end

%% find extraneous edges

% method: look for "triangles" of 3 voxels that are connected directly
% together. If the triangle is very flat (one edge is much longer than the
% sum of the other 2), the network would be just as efficient without that
% longer edge -- so trim it.


% loop through the upper triangle of the sparse network matrix
tic
for i = 1:nVox % loop through the rows    
    % optional shortcut for code dev
    if devFlag == 1 && length(trimCells) > 100
        break
    end
    
    if mod(i, 10) == 0
        fprintf('working through voxel %d / %d...\n', i, nVox)
    end
    
    for j = i+1:nVox % loop through the cols, as long as you stay in the upper triangle
        % optional shortcut for code dev
        if devFlag == 1 && length(trimCells) > 100
            break
        end
        
        % look for voxels that are directly connected to both i & j
        
        % multiply rows i and j -- voxels directly connected to both should have a 1 in the product:
        rowProd = sparseNet(i, :) .* sparseNet(j, :);
        
        % any candidate voxels?
        if ~isempty(nonzeros(rowProd)) 
            % get length of the edge betwen i & j:
            length_ij = voxDists(i, j);
            
            % get the voxels that are in a triangle with i & j, and only 
            % loop through those voxels, and focus on indices >= j+1 to 
            % avoid finding the same triangle multiple times:
            kVoxIdx = find(rowProd == 1);
            kVoxConsider = kVoxIdx(logical(kVoxIdx > j));

            % loop through the voxels (k) that are in a triangle with i & j
            for currk = 1:length(kVoxConsider)
                % optional shortcut for code dev
                if devFlag == 1 && length(trimCells) > 100
                    disp('breaking out into shorter version...')
                    break
                end
                
                k = kVoxConsider(currk);
                
                % get lengths of edges between i & k and j & k:
                length_ik = voxDists(i, k);
                length_jk = voxDists(j, k);
                
                % loop through the triangle's legs, and decide whether to 
                % trim any of them:
                legs = [length_ij, length_ik, length_jk];
                legIdx = [i, j; i, k; j, k]; % each row = edge
                legNames = {'ij', 'ik', 'jk'};
                trimLegs = [];
                for l = 1:length(legs)
                    % compare this leg to the other two
                    compIdx = ones(length(legs), 1);
                    compIdx(l) = 0;
                    
                    if sum(legs(logical(compIdx))) < allowed_inflation*legs(l)
                        trimLegs = [trimLegs, l];
                    end
                    
                end
                % if any of the legs seemed trimmable, delete the first
                % one arbitrarily:
                if ~isempty(trimLegs)
                    % which leg?
                    condemnedLeg = trimLegs(1);
                    trimVox1 = legIdx(condemnedLeg, 1);
                    trimVox2 = legIdx(condemnedLeg, 2);
                    
                    % store the location of this edge to trim it later:
                    rowMatch = find(prevRow == trimVox1);
                    colMatch = find(prevCol == trimVox2);
                    cellMatch = intersect(rowMatch, colMatch);
                    trimCells = [trimCells, cellMatch];
                    
                    % do it symmetrically:
                    rowMatchMirror = find(prevRow == trimVox2);
                    colMatchMirror = find(prevCol == trimVox1);
                    cellMatchMirror = intersect(rowMatchMirror, colMatchMirror);
                    trimCells = [trimCells, cellMatchMirror];                    
                    
                end
                
            end
            
        end
    end
end

disp('Done identifying edges to trim.')

%% make the final new sparseNet.

fprintf('Making the new network...\n')
% setup to re-build the updated sparse matrix (should be faster than
% individually setting values corresponding to deleted edges to 0)
newRow = prevRow;
newCol = prevCol;

% remove the extraneous edges:
newRow(trimCells) = [];
newCol(trimCells) = [];
newVals = ones(length(newRow), 1);

% make the final new sparseNet.
updatedNet = sparse(newRow, newCol, newVals, nVox, nVox);

% get rid of potential duplicates (probably unnecessary, but just to be safe):
updatedNet = spones(updatedNet);

%% update the weighted net

% make a new sparse net with the distances for the spared edges:
sparedEdges = find(updatedNet);
newWeights = voxDists(sparedEdges);
updatedWeightedNet = sparse(newRow, newCol, newWeights, nVox, nVox);

fprintf('...finished trimming edges in %d minutes!\n', toc/60)

end