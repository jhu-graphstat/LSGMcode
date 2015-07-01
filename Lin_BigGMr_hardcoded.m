function [match, clust_labels] = Lin_BigGMr_hardcoded(A, B, s, numdim, max_clust_size, max_seeds, diag_aug, proj_sphere, convex_init)
%input: A, B: adjacency matrices for the two graphs to match
%		s: the number of seeds
%        seeds are assumed to be vertices 1,2,3,...,s in both graphs
%		numdim: the embedding dimension
%		max_clust_size: the maximum number of vertices in a cluster
%           this is limited by the maximum sized graph that can be handled
%           via SGM (usually, this is ~1000-max_seeds)
%       max_seeds: the maximum number of seeds to use in matching the clusters.
%                   if SGM can handle graphs of size ~1000, then this would
%                   be ~1000-max_clust_size
%                   note that in the event of unbalanced clusters, it may
%                   be better to choose as many seeds as possible for each
%                   cluster.
%		diag_aug: if it is 1, we diagonally augment the matrices. if it is 0, we do not
%		proj_sphere: if it is 1, we project the embedding onto the sphere for subsequent clustering
%		convex_init: if it is 1, we start SGM with convex initialization 
%output: match: returns the permutation of the vertices, i.e., alignment of vertices.

% For example: BigGMr_hardcoded( AA,BB,m, numdim, max_clust,max_seeds, 0,0,0) is
% plain LSGM 
% BigGMr_hardcoded( AA,BB,m, numdim, max_clust,max_seeds, 1,1,1) is LSGM using diag. augmentation,
%       spherical projection, and convex_init. 
% Empirically, spherical projection improves accuracy significantly.

start = tic;
sumn = length(A)-s;
numclust = ceil(sumn/max_clust_size)+1;

% maxmium number of seeds to use
s_max = min(max_seeds,s);

% show output
show_output = true;

% choice for diagonal augmentation.
if diag_aug == 1
    % diagonal augmentation
	A = A + diag(sum(A,2)/(length(A)-1));
	B = B + diag(sum(B,2)/(length(B)-1));
else
    % no diagonal augmentation
	A = A;
	B = B;
end
% perform embedding
startt = tic;
embedMethod = 'adjacency'; % choose between adjacency and laplacian
if strcmp(embedMethod,'adjacency')
    [XA, XB] = spectralEmbed(A, B, numdim);
elseif strcmp(embedMethod,'laplacian')
    A = diag(1./sqrt(sum(A,2)))*A*diag(1./sqrt(sum(A,1)));
    B = diag(1./sqrt(sum(B,2)))*B*diag(1./sqrt(sum(B,1)));
    [XA, XB] = spectralEmbed(A, B, numdim);
else
    disp('Error: choose between adjacency and laplacian.')
end
if show_output
	fprintf( 'done projection: %f\n', toc(startt) );
end

% compute procrusties othogonal projection (on the seed vertices)
startt = tic;
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(sumn+s,1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
if show_output
	fprintf( 'done procrusties: %f\n', toc(startt) );
end

% choice for projection onto the sphere
if proj_sphere == 1
    % projection on the sphere
    XA = normr(XA);
	XB = normr(XB);
else
    % no projection
	XA = XA;
	XB = XB;
end

PD = pdist2(XA,XB);
Lin_DistMatching(PD);
[val,ind] = min(PD);
sum(ind == 1:size(XB,1))
%[Hmat,~] = YiCaoHungarian(-XA*XB');

figure(1)
scatter3(XA(s+1:end,1),XA(s+1:end,2),XA(s+1:end,3))
hold on
scatter3(XB(s+1:end,1),XB(s+1:end,2),XB(s+1:end,3))

% figure(1)
% scatter(XA(:,1),XA(:,2))
% hold on
% scatter(XB(:,1),XB(:,2))

% cluster (K-means) using the embedding
startt = tic;
XAXB=[XA;XB];
seedsA = 1:s;
seedsB = s+sumn+1:2*s+sumn;
nonseedsA = s+1:s+sumn;
nonseedsB = s+sumn+ s+1:2*(s+sumn);
nonseeds = [nonseedsA, nonseedsB];
[IDX, centroid, Dis] = kmeansAlgr(XAXB, numclust, nonseeds);
if show_output
	fprintf( 'done clustering: %f\n', toc(startt) );
end

figure(2)
ID_tmp = reshape(IDX,[length(IDX)/2,2]);
for i=1:numclust
    ai = XA(find(ID_tmp(:,1)==i),:);
    bi = XB(find(ID_tmp(:,2)==i),:);
    size(ai)
    size(bi)
    subplot(1,numclust,i)
    scatter3(ai(:,1),ai(:,2),ai(:,3))
    hold on
    scatter3(bi(:,1),bi(:,2),bi(:,3))
end

% fixing cluster sizes to be equal in both graphs
% we ensure the cluster sizes are consistent
% this can probably be done much more efficiently
startt = tic;
IDXs = IDX(seedsA);
IDXA = IDX(nonseedsA);
DisA = Dis(nonseedsA,:);
IDXB = IDX(nonseedsB);
DisB = Dis(nonseedsB,:);

clear IDX Dis

clustsizesA=zeros(numclust,1);
clustsizesB=zeros(numclust,1);
for i=1:numclust
   iiAi=find(IDXA==i);
   iiBi=find(IDXB==i);
   clustsizesA(i)=length(iiAi);
   clustsizesB(i)=length(iiBi);
end

clustsizes=round( (clustsizesA+clustsizesB)/2  );

% add or subtract 1 from largest cluster sizes
[val,ind] = sort(clustsizes, 1, 'descend');
temp = sumn -sum(clustsizes);
mask = ind(1:abs(temp));
clustsizes(mask) = clustsizes(mask) +sign(temp);

% store subgraphs induced by the clustering result
% this is used later for seeded graph matching within each cluster
pieceA_ = cell(numclust,1);
pieceB_ = cell(numclust,1);
gA_ = cell(numclust,1);
gB_ = cell(numclust,1);

disp('number of clusters')
numclust

for i=1:numclust
	if clustsizes(i) == 0
		continue
	end
	% sort by distance
	[~, TA] = sort(DisA(:,i));
	[~, TB] = sort(DisB(:,i));
	% use closest vertices
    gA = TA(1:clustsizes(i));
    gB = TB(1:clustsizes(i));
    % make these vertices unselectable for next time
	DisA(gA,:) = inf;
	DisB(gB,:) = inf;
    
    gA = gA +s;
    gB = gB +s;
	
    gAaug = [1:s gA'];
    gBaug = [1:s gB'];
    
    pieceA=A(gAaug,gAaug);
    pieceB=B(gBaug,gBaug);
    
    pieceA_{i} = pieceA;
    pieceB_{i} = pieceB;
    
    gA_{i} = gA;
    gB_{i} = gB;
end

clear A B DisA DisB

% perform graph matching in parallel
match = zeros(s+sumn,numclust);
clust_labels_ = zeros(s+sumn,2,numclust);
matching_init_time = toc(startt);
matching_time = 0;
parfor i = 1:numclust
	startt = tic;

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
	time_r = 0;
    % perform graph match
	if length(gA) > max_clust_size
        disp('lin!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		startr = tic;
        %'recurse'
		% cluster to large, match recursively
		[ord ii] = Lin_BigGMr_hardcoded( pieceA, pieceB, s, numdim, max_clust_size, max_seeds, diag_aug, proj_sphere, convex_init);
	    ii = ii(s+1:end,:);
	    
		ord = ord(s+1:end) -s;
		time_r = toc(startr);
		fprintf('recursion time: %f\n', time_r);
	else
  	% select seeds for SGM: active seed selection
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

			seeds(seed) = seedcans(ind);
			seedcans = seedcans([1:ind-1, ind+1:end]);
		end

		ind = [seeds, s+1:length(pieceA)];

		temp = zeros(s+sumn,2);
		
        % choice for initialization in graph matching
		if convex_init == 1
            % convex initialization
			ord = seedgraphmatchell2(pieceA(ind, ind), pieceB(ind, ind), s_max, 0);	
        else
            %barycenter
			ord = seedgraphmatchell2(pieceA(ind, ind), pieceB(ind, ind), s_max); 
		end
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
    matching_time = toc(startt)-time_r + matching_time;
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
	fprintf( 'done matching: %f\n', matching_init_time+matching_time);
	fprintf( 'total time: %f\n', toc(start) );
end

end

