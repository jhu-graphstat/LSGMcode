clear all
% load the data
name = 'w5000';
load([name,'.mat'])
A = spfun(@logWeight,A);
B = spfun(@logWeight,B);
num_ht = double(num_ht);
num_at = double(num_at);

% initialize parameters
N = length(A);  % number of vertices
num_seeds = num_ht+num_at;  % number of seeds 
max_clust = 500;   % max cluster size
numdim = 20;    % maximum dimension considered for the embedding

% one-hop neighbors
seedinds = 1:num_seeds;
[pA,pB] = seedCoverage(seedinds,A,B);
fprintf( 'Seed coverage (graph A): %f\n', pA )
fprintf( 'Seed coverage (graph B): %f\n', pB )

%% compute accuracy a
[match,clust_labels_,num_partition] = Lin_BigGMr_hashtags(A,B,num_ht,num_seeds, numdim, max_clust, 0,0,0);
% Aperm_inv = zeros(1,N);
% Aperm_inv(Aperm) = 1:N;
% Bperm_inv = zeros(1,N);
% Bperm_inv(Bperm) = 1:N;
% num_clust = size(clust_labels_,3);
% results_ = zeros(N,num_clust);
% results = zeros(1,num_clust);
% for i = 1:num_clust
%     clust_labels = clust_labels_(:,:,i);
%     %re-arrange clust_labels
%     results_(find(clust_labels(Aperm_inv,1)+clust_labels(Bperm_inv,2) == 2),i) = 1;
%     results(i) = sum(clust_labels(Aperm_inv,1)+clust_labels(Bperm_inv,2) == 2)/min([sum(clust_labels(:,1)),sum(clust_labels(:,2))]);
% end
% fprintf('number of nodes in each partition')
num_partition
% results
% fprintf('oracle accuracy %f\n',sum(sum(results_,2)>=1)/(N-num_seeds))

% % computing matching accuracy
% B_match = zeros(num_clust,N-num_seeds);
% results = zeros(1,num_clust);
% for i = 1:num_clust
%     ind = find(match(:,i));
%     b = Bperm(match(ind,i));
%     B_match(i,ind-num_seeds) = b; % zero for un-matched nodes
%     match_ab = (Aperm(ind)==b);
%     results(i) = mean(match_ab);
% end
% results
% tmp = [Aperm(num_seeds+1:N);sum(B_match)];
% mean(tmp(1,:)==tmp(2,:))


%% graph matching on the entire graph
% convex_init = 1;
% s_max = round((size(AA,1)-num_ht));
% seeds = selectSeeds(AA,BB,num_seeds,s_max);
% ind = [seeds, num_seeds+1:length(AA)]; % indices for both seeds and nodes to be matched
% % choice for initialization in graph matching
% if convex_init == 1
%     % convex initialization
%     ord = seedgraphmatchell2(AA(ind, ind), BB(ind, ind), s_max, 0);	
% else
%     %barycenter
%     ord = seedgraphmatchell2(AA(ind, ind), BB(ind, ind), s_max); 
% end
% ord = ord(s_max+1:end)-s_max;
%     
% % save results
% gB = num_seeds+1:length(A);
% mean(Aperm(num_seeds+1:end) == Bperm(gB(ord)))
