function updatedNet = trimEdges(sparseNet, voxDists, allowed_inflation)
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

fprintf('\ntrimming edges...\n')

% how many voxels?
nVox = size(sparseNet, 1);

% setup to store the edges being trimmed (should be faster than updating
% the sparse net directly, even using some tricks that still require some
% really long vectors):
trimRow = [];
trimCol = [];

%% find extraneous edges

% method: look for "triangles" of 3 voxels that are connected directly
% together. If the triangle is very flat (one edge is much longer than the
% sum of the other 2), the network would be just as efficient without that
% longer edge -- so trim it.


% loop through the upper triangle of the sparse network matrix
tic
for i = 1:nVox % loop through the rows
    if mod(i, 100) == 0
        fprintf('working through voxel %d / %d...\n', i, nVox)
    end
    
    for j = i+1:nVox % loop through the cols, as long as you stay in the upper triangle
        % look for voxels that are directly connected to both i & j
        
        % multiply rows i and j -- voxels directly connected to both should have a 1 in the product:
        rowProd = sparseNet(i, :) .* sparseNet(j, :);
        
        % any candidate voxels?
        if ~isempty(nonzeros(rowProd)) 
            % get length of the edge betwen i & j:
            length_ij = voxDists(i, j);
            
            % focus on indices >= j + 1 to avoid finding the same triangle
            % multiple times
            for k = j+1:nVox
                if rowProd(k) == 1 % voxel k is in a triangle with j & k
                    
                    % get lengths of edges between i & k and j & k:
                    length_ik = voxDists(i, k);
                    length_jk = voxDists(j, k);
                    
                    % loop through the legs, and decide whether to trim any of
                    % them:
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
                        trimRow = [trimRow, trimVox1]; % don't need to make this symmetrical yet.
                        trimCol = [trimCol, trimVox2];

                    end
                end
            end
            
        end
    end
end

fprintf('...finished trimming edges in %d sec!\n', toc)

%% make the final new sparseNet.

keyboard

% setup to re-build the updated sparse matrix (should be faster than
% individually setting values corresponding to deleted edges to 0)
[rowIdx, colIdx, vals] = find(sparseNet); % pull out the value & location of all non-zero elements


% [] trim the trimmed edges (*symmetrically*)
% [] make a new net

% [] check this
updatedNet = sparse(edgeRow, edgeCol, edgeVals, nVox, nVox);


end