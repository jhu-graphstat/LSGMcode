function [ match, clust_labels ] = BigGM( A, B, s, numdim, numclust, embedAlg, clustAlg, graphMatchAlg)
% BigGM Large Seeded Graph Matching function
% [ match, clust_labels ]  = BigGM( A, B, s, numdim, numclust, embedAlg, clustAlg, graphMatchAlg)
% A generic function to perform graph matching on large graphs via a
% divided and conquer strategy. The steps of the procedure are as follows:
% (1) each graph is embedded into Euclidean space,
% (2) the embeddings are aligned using the seed vertices,
% (3) the points/vertices are clustered 

switch nargin
    case 3
        warning(['Embedding dimension dimension and number of clusters '...
            'not supplied; both set to a value of 6. This is silly.'])
        numdim = 6; % A silly default
        numclust = 6; % Equally silly default
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeans0; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite 
    case 4
        warning(['Number of clusters '...
            'not supplied; default value is 6. This is silly.'])
        numclust = 6; % Equally silly default
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeans0; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
    case 5
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmean0; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
    case 6
        clustAlg = @kmeans0; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
    case 7
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
end

start = tic;
sumn = length(A)-s;

%% perform embedding
startt = tic;
[XA XB] = embedAlg(A, B, numdim);
fprintf( 'done projection: %f\n', toc(startt) );

%% compute procrusties othogonal projection (on the seed vertices)
startt = tic;
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(sumn+s,1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
fprintf( 'done procrusties: %f\n', toc(startt) );


%% cluster using the embedding
startt = tic;
XAXB=[XA;XB];
nonseedsA = s+1:s+sumn;
nonseedsB = s+sumn+ s+1:2*(s+sumn);
[IDX, centroid, Dis] = clustAlg(XAXB, numclust);
fprintf( 'done clustering: %f\n', toc(startt) );

%% fix cluster sizes to be equal in both graphs

[pieceA_,pieceB_,gA_,gB_] = fixClusterSize(A,B,IDX, Dis, numclust,nonseedsA, nonseedsB);
clear IDX Dis

%% perform graph matching in parallel
match = zeros(s+sumn,numclust);
parfor i = 1:numclust
	% load subgraph adjacency matrix
	gA = gA_{i};
	gB = gB_{i};
	pieceA = pieceA_{i};
	pieceB = pieceB_{i};
	
	% if graph is empty
	if length(gA)==0
		continue;
	end
	
    % perform graph match
	ord = graphMatchAlg(pieceA, pieceB, s);
	ord = ord(s+1:end)-s
	
	% save results
	temp = zeros(s+sumn,1);
	temp(gA) = gB(ord);
    match(:,i) = temp;
end
% combine results
match(1:s,1) = 1:s;
match = sum(match,2)';

fprintf( 'done matching: %f\n', toc(startt) );
fprintf( 'total time: %f\n', toc(start) );

end
