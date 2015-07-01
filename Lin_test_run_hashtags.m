clear all

%name = 'graphs525';
%num_at = 525;

name = 'graphs4181';
num_at = 4181;

load([name,'.mat'])
A = spfun(@logWeight,A);
B = spfun(@logWeight,B);

%% set parameters
% number of vertices
N = length(A);
num_ht = N-num_at;

% number of seeds
num_seeds = num_ht+500;  % need agent seeds graph partition step
 
% max cluster size
max_clust = 800;

% maximum dimension considered for the embedding
numdim = 20;

%% random permutating the nodes
tmp = num_ht+1:length(A);
Aperm = [1:num_ht,tmp(randperm(num_at))]; % Hashtags are seeds and it is ordered.
Bperm = [1:num_ht,tmp(randperm(num_at))];

% select seeds
seedinds = Aperm(1:num_seeds);
[pA,pB] = seedCoverage(seedinds,A,B);
fprintf( 'Seed coverage (graph A): %f\n', pA )
fprintf( 'Seed coverage (graph B): %f\n', pB )

% move vertices in B to match seedinds
for seed_ind = num_ht+1:num_seeds
	seed = seedinds(seed_ind);
	Bperm(Bperm==seed) = Bperm(seed_ind);
	Bperm(seed_ind) = seed;
end

%Bend = Bperm(num_seeds+1:end);
AA = A(Aperm,Aperm);
BB = B(Bperm,Bperm);
% 
%% compute accuracy a
strategy = 1;  % 1: recursive kmean, 2: infomap community detection (perform kmean once whenever the subgraph graph has only one community) 3: alternating between infomap and kmean

[match,clust_labels_,num_partition] = Lin_BigGMr_hashtags(AA,BB,num_ht,num_seeds, numdim, max_clust,strategy,0,0,0);
Aperm_inv = zeros(1,N);
Aperm_inv(Aperm) = 1:N;
Bperm_inv = zeros(1,N);
Bperm_inv(Bperm) = 1:N;
num_clust = size(match,2);
results_ = zeros(N,num_clust);
results = zeros(1,num_clust);
for i = 1:num_clust
    %re-arrange clust_labels
    results_(find(clust_labels_{1}(Aperm_inv,i)+clust_labels_{2}(Bperm_inv,i) == 2),i) = 1;
    results(i) = sum(clust_labels_{1}(Aperm_inv,i)+clust_labels_{2}(Bperm_inv,i) == 2)/min([sum(clust_labels_{1}(:,i)),sum(clust_labels_{2}(:,i))]);
end
fprintf('number of nodes in each partition')
num_partition
results
fprintf('oracle accuracy %f\n',sum(sum(results_,2)>=1)/(N-num_seeds))

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

