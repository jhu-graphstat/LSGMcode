function [ match, clust_labels ] = BigGMr( A, B, s, numdim, max_clust_size, embedAlg, clustAlg, graphMatchAlg)
% BigGM Large Seeded Graph Matching function
% [ match, clust_labels ]  = BigGM( A, B, s, numdim, numclust, embedAlg, clustAlg, graphMatchAlg)
% A generic function to perform graph matching on large graphs via a
% divided and conquer strategy. The steps of the procedure are as follows:
% (1) each graph is embedded into Euclidean space,
% (2) the embeddings are aligned using the seed vertices,
% (3) the points/vertices are clustered 
% (4) in parallel, match subgraphs corresponding to each cluster
%       if the cluster size is to big the algorithm recurses to 
%       step 1 on the subgraphs

switch nargin
    case 3
        warning(['Embedding dimension dimension and number of clusters '...
            'not supplied; both set to a value of 6. This is silly.'])
        numdim = 6; % A silly default
        numclust = 6; % Equally silly default
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite 
    case 4
        warning(['Number of clusters '...
            'not supplied; default value is 6. This is silly.'])
        numclust = 6; % Equally silly default
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
    case 5
        embedAlg = @spectralEmbed; % Use the adjacency spectral embedding
        clustAlg = @kmeanAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
    case 6
        clustAlg = @kmeansAlg; % Use kmeans
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
    case 7
        graphMatchAlg = @seedgraphmatchell2; % Our favorite
end

start = tic;
sumn = length(A)-s;

% number of clusters determined by max cluster size
numclust = ceil(sumn/max_clust_size);
% maxmium number of seeds to use
s_max = min(200,s);
% show output
show_output = false;

%% perform embedding
startt = tic;
[XA XB] = embedAlg(A, B, numdim);
if show_output
	fprintf( 'done projection: %f\n', toc(startt) );
end

%% compute procrusties othogonal projection (on the seed vertices)
startt = tic;
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(sumn+s,1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
if show_output
	fprintf( 'done procrusties: %f\n', toc(startt) );
end

%% cluster using the embedding
startt = tic;
XAXB=[XA;XB];
nonseedsA = s+1:s+sumn;
nonseedsB = s+sumn+ s+1:2*(s+sumn);
nonseeds = [nonseedsA, nonseedsB];
[IDX, centroid, Dis] = clustAlg(XAXB, numclust);
if show_output
	fprintf( 'done clustering: %f\n', toc(startt) );
end
%% fix cluster sizes to be equal in both graphs

[pieceA_,pieceB_,gA_,gB_] = fixClusterSize(A,B,IDX, Dis, numclust,nonseedsA, nonseedsB);
clear IDX Dis

%% perform graph matching in parallel
match = zeros(s+sumn,numclust);
clust_labels_ = zeros(s+sumn,2,numclust);

% Need a for loop here because of recursion, :(
for i = 1:numclust
	
	ii = [1,1];
	
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
	if length(gA) > max_clust_size
		startr = tic;
%		'recurse'
		% cluster to large, match recursively
		[ord ii] = BigGMr( pieceA, pieceB, s, numdim, max_clust_size, embedAlg, clustAlg, graphMatchAlg);
	    ii = ii(s+1:end,:);
	    
		ord = ord(s+1:end) -s;
		time_r = toc(startr);
		fprintf('recursion time: %f\n', time_r);
	else
%		'cluster'

    	% select seeds for SGM
    	seedsA = pieceA(1:s, s+1:end);
    	seedsB = pieceB(1:s, s+1:end);
    	
    	avgdegA = mean(seedsA,2);
    	avgdegB = mean(seedsB,2);
    	seedcans = 1:s;
    	% compute first seed
		[val, ind] = min(abs(avgdegA(seedcans) -.5) +...
						 abs(avgdegB(seedcans) -.5));
    	
    	seeds = zeros(1,s_max);
    	seeds(1) = ind;
    	
    	% set of remaining seed candidates
		seedcans = setdiff(seedcans, ind);
%    	pA = seedsA(ind,:)/2;
%    	pB = seedsB(ind,:)/2;
		p = [seedsA(ind,:), seedsB(ind,:)]/2;
    	for seed = 2:s_max
    		temp = repmat(p, [length(seedcans),1]);
    		p_ = [seedsA(seedcans,:), seedsB(seedcans,:)]/seed+temp;
    		p_1 = 1-p_;
    		temp = p_.*log(p_) +p_1.*log(p_1);
    		temp(find(isnan(temp))) = 0;
    		H = -sum(temp,2);
    	
    		[val, ind] = max(H);
    		p = p*seed/(seed+1);
    		
%			[val, ind] = min(abs(avgdegA(seedcans) -.5) +...
%							 abs(avgdegB(seedcans) -.5));
						 
%			[val, ind] = min(abs(avgdegA(seedcans) -.5 +penaltyA) +...
%							 abs(avgdegB(seedcans) -.5 +penaltyB));
			seeds(seed) = seedcans(ind);
			seedcans = seedcans([1:ind-1, ind+1:end]);
		end
%		seeds = ind(1:num_seeds);
		ind = [seeds, s+1:length(pieceA)];

	
%		% select seeds from the same cluster
%		%  and randomly select remaining seeds from other clusters
%		inds = zeros(s_max+size(pieceA,1)-s,1);
%		n1 = min(sum(IDXs==i),s_max);
%		inds(1:s_max) = [randsample( find(IDXs==i)', n1 ) , randsample(find(IDXs~=i)', s_max-n1) ];
%		% add remaining indices
%		inds(s_max+1:end) = [ s+1:size(pieceA,1)];
		temp = zeros(s+sumn,2);
		
		ord = graphMatchAlg(pieceA(ind, ind), pieceB(ind, ind), s_max);
		ord = ord(s_max+1:end)-s_max;
	end
	
	% save cluster labels
	temp = zeros(s+sumn,2);
	temp(gA, 1) = ii(:,1);
	temp(gB, 2) = ii(:,2);
	clust_labels_(:,:,i) = temp;
	
	% save results
	temp = zeros(s+sumn,1);
	temp(gA) = gB(ord);
    match(:,i) = temp;
    
end
% combine results
match(1:s,1) = 1:s;
match = sum(match,2)';

clust_labels = clust_labels_(:,:,1);
maxx = max(max(clust_labels));
for i = 2:numclust
	temp = clust_labels_(:,:,i);
	mask = temp ~= 0;
	clust_labels(mask) = clust_labels(mask) +temp(mask) +maxx;
	maxx = max(max(clust_labels));
end

if show_output
	fprintf( 'done matching: %f\n', toc(startt)-time_r );
	fprintf( 'total time: %f\n', toc(start) );
end

end

