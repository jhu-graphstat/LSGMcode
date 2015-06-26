function subA = prune(A,seeds)
% Looks at the connected components of the graph given by the adjacency
% matrix A, and returns a cell array of adjacency matrices corresponding to
% the connected components that contain a seed.

% REQUIRES A FUNCTION TO COMPUTE THE CONNECTED COMPONENTS OF A GRAPH WITH
% THE NAME components

[ci, ~] = components(A);
seedComponent = unique(ci(seeds));
n = length(seedComponent);

subA = cell(1,n);
for iComponent = 1:n
    subA{iComponent} = find(ci == seedComponent(iComponent));
end

end