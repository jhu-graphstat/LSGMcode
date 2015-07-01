function [match,clust_labels_,num_partition] = Lin_BigGMr_hashtags(A, B, num_ht, s, numdim, max_clust_size, strategy, diag_aug, proj_sphere, convex_init)
%%input: A, B: adjacency matrices for the two graphs to match
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

%%
% choice for diagonal augmentation.
if diag_aug == 1
    % diagonal augmentation
	A = A + diag(0.5*(sum(A,1)+sum(A,2)')/(length(A)-1));
	B = B + diag(0.5*(sum(B,1)+sum(B,2)')/(length(B)-1));
else
    % no diagonal augmentation
	A = A;
	B = B;
end

ind_start = num_ht;
num_seeds = s-ind_start; % number of non-hashtag seeds
AA = A(ind_start+1:end,ind_start+1:end);
BB = B(ind_start+1:end,ind_start+1:end);

if strategy == 1
    [gA_,gB_] = embeddingKmean(AA,BB,num_seeds,max_clust_size,numdim);
elseif strategy == 2
    [gA_,gB_] = infomapClust(AA,BB,num_seeds,max_clust_size,numdim,2,1);
else
    [gA_,gB_] = infomapKmean(AA,BB,num_seeds,max_clust_size,numdim,2,1);
end
numclust = length(gA_);
num_partition = [cellfun(@length,gA_) ; cellfun(@length,gB_)];

% check the oracle accuracy
clust_labels_ = {};
clust_labels_{1} = zeros(length(A),numclust);
clust_labels_{2} = zeros(length(B),numclust);
for i = 1:numclust
    % save cluster labels
    tempA = zeros(length(A),1);
    tempA(gA_{i}+ind_start) = 1;
    tempB = zeros(length(B),1);
    tempB(gB_{i}+ind_start) = 1;
    clust_labels_{1}(:,i) = tempA; % indicator function
    clust_labels_{2}(:,i) = tempB;
end

%% perform graph matching in parallel
match = zeros(length(A),numclust); % list of nodes in B matched to nodes in A 
% for i = 1:numclust	% can  be replaced with parfor along with some minor changes
%     fprintf('cluster number %f',i)
% 	% load subgraph adjacency matrix
% 	gA = gA_{i}+ind_start;
% 	gB = gB_{i}+ind_start;
%     
%     gAaug = [1:s, gA];
%     gBaug = [1:s, gB];
%     pieceA=A(gAaug,gAaug);
%     pieceB=B(gBaug,gBaug);
%     
% 	% if graph is empty
%     if isempty(gA)
% 		continue;
%     end
%     % perform graph match
%     s_max = max(round((size(pieceA,1)-s)*2),50);
%     fprintf('number of s_max %d\n',s_max)
%     seeds = selectSeeds(pieceA,pieceB,s,s_max);
%     seeds_ht = sum(seeds<=num_ht);
%     seeds_at = sum(seeds>num_ht);
%     fprintf('number of hashtag seeds %d\n', seeds_ht)
%     fprintf('number of agent seeds %d\n', seeds_at)
% 	indA = [seeds, s+1:length(pieceA)]; % indices for both seeds and nodes to be matched
%     indB = [seeds, s+1:length(pieceB)];
%     AA = pieceA(indA,indA);
%     BB = pieceB(indB,indB);
%     % check the number of nodes from each graph
%     len_A = length(AA);
%     len_B = length(BB);
%     % if the number of nodes are different, place dummy nodes in the
%     % smaller graph
%     if len_A > len_B
%         BB = blkdiag(spfun(@addOne,BB)-ones(len_B), zeros(len_A-len_B));
%         AA = spfun(@addOne,AA)-ones(len_A);
%         %BB = blkdiag(BB, zeros(len_A-len_B));    
%     elseif len_B > len_A
%         AA = blkdiag(spfun(@addOne,AA)-ones(len_A), zeros(len_B-len_A));
%         BB = spfun(@addOne,BB)-ones(len_B);
%         %AA = blkdiag(AA, zeros(len_B-len_A));
%     end
%     % SGM alg.
% 	ord = seedgraphmatchell2(AA, BB, s_max,convex_init); 
% 	ord = ord(s_max+1:end)-s_max;
% 
%     % save results
%     temp = zeros(length(A),1);
%     if len_A > len_B
%         ind = find(ord<=len_B-s_max);
%         temp(gA(ind)) = gB(ord(ind));
%     elseif len_B > len_A
%         temp(gA) = gB(ord(1:len_A-s_max));
%     else
%         temp(gA) = gB(ord);
%     end
%     match(:,i) = temp;
% end

function b = addOne(a)
b = a+1;