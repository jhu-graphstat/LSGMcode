function [pieceA_,pieceB_,gA_,gB_] = processClusters(A,B,IDX, nonseedsA, nonseedsB)
% [pieceA_,pieceB_,gA_,gB_] = processClusters(A,B,IDX, numclust,nonseedsA, nonseedsB)
% Returns four cell arrays with (respectively) the adjacency submatrices
% for graphs A and B corresponding to each cluster, and arrays of indices
% of the vertices in graphs A and B in each cluster.
%
% INPUTS:       A, B : (sparse matrix) adjacency matrices
%                IDX : array with cluster assignment for all vertices in
%                      both graphs
%          nonseedsA : array of indices of nonseed vertices from graph A
%          nonseedsB : array of indices of nonseed vertices from graph B
% OUTPUTS:   pieceA_ : cell array whose i-th entry contains the adjacency
%                      submatrix for vertices in cluster i and graph A
%            pieceB_ : cell array whose i-th entry contains the adjacency
%                      submatrix for vertices in cluster i and graph B
%                gA_ : cell array whose i-th entry contains an array of the
%                      indices for vertices in the i-th cluster and graph A
%                gB_ : cell array whose i-th entry contains an array of the
%                      indices for vertices in the i-th cluster and graph B

nANonseeds = length(nonseedsA);
s = length(A)-nANonseeds;

numclust = max(IDX);


IDXA = IDX(nonseedsA);
IDXB = IDX(nonseedsB);

% save clusters(subgraphs) to match (so it can be parallelized)
pieceA_ = cell(numclust,1);
pieceB_ = cell(numclust,1);
gA_ = cell(numclust,1);
gB_ = cell(numclust,1);

for i=1:numclust
    % find i-th clusters in both graphs
    iACluster=find(IDXA==i);
    iBCluster=find(IDXB==i);
	if isempty(iACluster) || isempty(iBCluster)
		continue
	end
   
    gA = iACluster +s;
    gB = iBCluster +s;
	
    gAaug = [1:s gA'];
    gBaug = [1:s gB'];
    
    % Get adjacency submatrix
    pieceA=A(gAaug,gAaug);
    pieceB=B(gBaug,gBaug);
    
    nPieceA = length(pieceA);
    nPieceB = length(pieceB);
    
    if (nPieceA > 10000 || nPieceB > 10000)
        error('A cluster contains more than 10000 nodes: exiting to prevent memory issues')
    elseif (nPieceA > 1000 || nPieceB > 1000)
        warning('A cluster contains more than 100- nodes: expect performance decrease')
    end

    
    pieceA_{i} = pieceA;
    pieceB_{i} = pieceB;
    
    gA_{i} = gA;
    gB_{i} = gB;
end
