function [ match, clust_labels ] = BigGMr( A, B, s, numdim, max_clust_size, embedAlg, clustAlg, graphMatchAlg, topK)
% BigGM Large Seeded Graph Matching function
% [ match, clust_labels ] = BigGMr( A, B, s, numdim, max_clust_size, embedAlg, clustAlg, graphMatchAlg, topK)
% A generic function to perform graph matching on large graphs via a
% divided and conquer strategy. The steps of the procedure are as follows:
% (1) each graph is embedded into Euclidean space,
% (2) the embeddings are aligned using the seed vertices,
% (3) the points/vertices are clustered 
% (4) in parallel, match subgraphs corresponding to each cluster
%       if the cluster size is to big the algorithm recurses to 
%       step 1 on the subgraphs
%
% INPUTS:          A, B : adjacency matrices for each graph
%                     s : number of seeds
%                numdim : dimension to be used in the embedding algoritm
%        max_clust_size : maximum number of points in each cluster allowed
%                         when clustering the embedded points. This is 
%                         used to calculate the number of clusters used.
%              embedAlg : algorithm used to embed the adjacency matrices
%              clustAlg : algorithm used to cluster the embedded adjacency
%                         matrices
%         graphMatchAlg : algorithm used to match vertices from each graph
%                         in a cluster
%                  topK : boolean value set to TRUE when a list of matches
%                         is desired


switch nargin
    case 3
        warning(['Embedding dimension dimension and number of clusters '...
            'not supplied; both set to a value of 6. This is silly.'])
        numdim = 6; % A silly default
        numclust = 6; % Equally silly default
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite 
        topK = false;
    case 4
        warning(['Number of clusters '...
            'not supplied; default value is 6. This is silly.'])
        numclust = 6; % Equally silly default
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
        topK = false;
    case 5
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
        topK = false;
    case 6
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
        topK = false;
    case 7
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
        topK = false;
    case 8
        topK = false;
end

start = tic;
nANonseeds = length(A)-s;
nBNonseeds = length(B) - s;

sumn = max(nANonseeds, nBNonseeds);

% number of clusters determined by max cluster size
numclust = ceil(sumn/max_clust_size);
% maxmium number of seeds to use
s_max = min(200,s);
% show output
 show_output = true;

% Log transform count data

log1 = @(x) log(x+1);

A = spfun(log1, A);
B = spfun(log1, B);

% Use regularized graph laplacian (Qin and Rohe 2012)

LA = regularizedLaplacian(A);
LB = regularizedLaplacian(B);

[IDX, nonseedsA, nonseedsB] = jointCluster(LA, LB, s, numdim, max_clust_size);


%% fix cluster sizes to be equal in both graphs
if topK == true
    [pieceA_,pieceB_,gA_,gB_] = processClusters(A,B,IDX, nonseedsA, nonseedsB);
else
    [pieceA_,pieceB_,gA_,gB_] = fixClusterSize(A,B,IDX, Dis, numclust,nonseedsA, nonseedsB);
end
clear IDX Dis

%% perform graph matching in parallel
startt = tic;
if (topK == true)
    match = sparse(nAnonseeds + s, nBnonseeds + s); % entry i,j is probability vertex i in A is matched to vertex j in B
else
    match = zeros(s+sumn,numclust);
end
clust_labels_ = zeros(s+sumn,2,numclust);

% Need a for loop here because of recursion, :(
for i = 1:numclust
	
	ii = [1,1];
	
	% load subgraph adjacency matrix
	gA = gA_{i};
	gB = gB_{i};
	

	% if either graph is empty 
	if (isempty(gA) || isempty(gB))
		continue;
    end
    
    [pieceA, pieceB] = addDummyNodes(pieceA_{i}, pieceB_{i});
    
    ngA = length(gA); % gA and gB can be different sizes!
    ngB = length(gB);

%		'cluster'

    % select seeds for SGM
    seedsA = pieceA(1:s, s+1:end);
    seedsB = pieceB(1:s, s+1:end);
    seeds = activeSeedSelection(seedsA, seedsB, s_max);
   
    ind = [seeds, s+1:length(pieceA)];

    ord = graphMatchAlg(pieceA(ind, ind), pieceB(ind, ind), s_max, topK);
    if (topK == false)
        ord = ord(s_max+1:end)-s_max;
    end
	
	% save cluster labels
	temp = zeros(s+sumn,2);
    % Hack to handle different sized graphs
    if (length(ii) > 2)
        nA = length(gA);
        nB = length(gB);
        temp(gA, 1) = ii(1:nA,1);
        temp(gB, 2) = ii(1:nB,2);
    else
        temp(gA, 1) = ii(:,1);
        temp(gB, 2) = ii(:,2);
    end
	clust_labels_(:,:,i) = temp;
	

	% save results
    if (topK == true)

        %nSeeds = length(seeds);
        % A not very clever way to record the results

        % A simple way to record the Nonseeds
        for j = 1:ngA
            for k = 1:ngB
                match(gA(j), gB(k)) = ord(s_max + j, s_max + k);
            end
        end
    else
        temp = zeros(s+sumn,1);
        temp(gA) = gB(ord);
        match(:,i) = temp;
    end
    
    
end
% combine results
if (topK == true)
    for i = 1:s
        match(i,i) = 1; % seeds are always matched to themselves
    end
elseif (topK == false)
    match(1:s,1) = 1:s;
    match = sum(match,2)';
end

clust_labels = clust_labels_(:,:,1);
maxx = max(max(clust_labels));
for i = 2:numclust
	temp = clust_labels_(:,:,i);
	mask = temp ~= 0;
	clust_labels(mask) = clust_labels(mask) +temp(mask) +maxx;
	maxx = max(max(clust_labels));
end

if show_output
	fprintf( 'done matching: %f\n', toc(startt) );
	fprintf( 'total time: %f\n', toc(start) );
end


end

