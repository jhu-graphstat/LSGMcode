function [best_v, vertexObj] = getBestMatches(A,B,coor, seeds,which)
% Returns the vertices that are best matched by the given correspondence
% 
% best_v = getBestMatches(a,b,coor, seeds)
% 
% best_v is the nonSeeds sorted according to how much they contribute to
%   the objective function
%
% A,B are two adjancency matrices 
% coor is the predicted correspondence between vertices
% seeds is the list of seed vertices

    if nargin<5
        which = 'normalized';
    end


    % How many nonseeds are there
    numNodes = numel(coor);
    numNonSeeds = numNodes-numel(coor(seeds));
    
    A = A(coor,coor);

    % Get the vertex-wise objective function normalized
    switch(which)
        case 'normalized'
            vertexObj = (ones(1,numNodes)*(B-A).^2)./...
                (ones(1,numNodes)*(B+A)+1);
        case 'raw'
            vertexObj = (ones(1,numNodes)*(B-A).^2);
        case 'corr'
            vertexObj = arrayfun(...
                @(v) 1-nzcorr(full(B(:,v)),full(A(:,v))),...
                1:numNodes);
    end
    
     % Set seeds objective function to Inf
     vertexObj(seeds) = Inf;
     
     % sort the objective function
     [~, idx] = sort(vertexObj);
    
    % return the list of idx
    best_v = idx(~ismember(idx,seeds));
end

function c = nzcorr(a,b)
nz = (a>0) | (b>0);
if sum(nz) > 1
    c = corr(a(nz),b(nz));
else
    c = 1;
end
end