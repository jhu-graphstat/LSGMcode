
%matlabpool(12); % depends on how many cores you have on the machine. max = 12.

% Load graphes (adjacency matrices A and B)
%load('KKI-08-29.mat'); 

load('A_binary.mat')
load('B_binary.mat')
%A = 0.5*A+0.5*A';
%B = 0.5*B+0.5*B';

% name of save file
save_file_name = 'lsgm_data_lin_test_run.mat';

% max cluster size
max_clust = 300;

% number of seeds
num_seeds = 200;
max_seeds = 100; % the maximum number of seeds to use in matching the clusters

% number of vertices
N = size(A,1);

% maximum dimension considered for the embedding
numdim = 3;

% make experiments reproducable
%rng(0);
	
% index vector for nonseed vertices
m = num_seeds;
nonseeds = m+1:N;

% random permutating the nodes
Aperm = randperm(N);
Bperm = randperm(N);

% select seeds
seedingMethod = 'random'; % choose from one of the five choices: descend, ascend, middle, random,combined
a = full(sum(A,1));
[sortedValues,sortIndex] = sort(a(:),'descend');
if strcmp(seedingMethod,'random')
    seedinds = Aperm(1:num_seeds);
elseif strcmp(seedingMethod,'descend')
    seedinds = sortIndex(1:num_seeds); 
elseif strcmp(seedingMethod, 'ascend')
    seedinds = sortIndex(end-num_seeds+1:end); 
elseif strcmp(seedingMethod,'middle')
    n = floor(length(sortIndex)/2);
    seedinds = sortIndex(n-num_seeds/2+1:n+num_seeds/2);
elseif strcmp(seedingMethod,'combined')
    n = floor(length(sortIndex)/2);
    ind = [1:num_seeds/4,n-num_seeds/4+1:n+num_seeds/4,N-num_seeds/4+1:N];
    tmp = setdiff(Aperm,ind);
    ind = [ind,tmp(1:num_seeds/4)];
    seedinds = sortIndex(ind);
else
    disp('Error: seeding method is not correct.')
end
[pA,pB] = seedCoverage(seedinds,A,B);
fprintf( 'Seed coverage (graph A): %f\n', pA )
fprintf( 'Seed coverage (graph B): %f\n', pB )

% move vertices in A to match seedinds 
for seed_ind = 1:num_seeds
	seed = seedinds(seed_ind);
	Aperm(Aperm==seed) = Aperm(seed_ind);
	Aperm(seed_ind) = seed;
end
% move vertices in B to match seedinds
for seed_ind = 1:num_seeds
	seed = seedinds(seed_ind);
	Bperm(Bperm==seed) = Bperm(seed_ind);
	Bperm(seed_ind) = seed;
end
%Bend = Bperm(num_seeds+1:end);
AA = A(Aperm,Aperm);
BB = B(Bperm,Bperm);
    
% lsgm
%rng(0);

[match, clust_labels] = Lin_BigGMr_hardcoded( AA,BB,m, numdim, max_clust, max_seeds, 1,1,1);
match_ab = (Aperm(nonseeds)==Bperm(match(nonseeds)));
% check the degree of matched nodes
Aperm_ns = Aperm(nonseeds);
a = Aperm_ns(find(match_ab));
b = Aperm_ns(find(~match_ab));
% figure(3)
% hold on
% plot(sort(sum(A(a,:),2)))
% plot(sort(sum(A(b,:),2)))

mean(match_ab)

% oracle accuracy
Aperm_inv = 1:N;
Aperm_inv(Aperm) = 1:N;
Bperm_inv = 1:N;
Bperm_inv(Bperm) = 1:N;
[max_clust mean(clust_labels(Aperm_inv(nonseeds),1)==clust_labels(Bperm_inv(nonseeds),2))]

%save(save_file_name, 'num_seeds_vec', 'num_max_seeds', 'numdim', 'max_clust_vec', 'acc', 'runtime');
