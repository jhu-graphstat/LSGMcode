%matlabpool(12); % depends on how many cores you have on the machine. max = 12.

% Load graphes (adjacency matrices A and B)
load('KKI-08-29.mat'); 

% name of save file
save_file_name = 'lsgm_data_08-29_test_run.mat';

% simulation parameters
num_runs = 5;
num_exp = 2;

% max cluster size
max_clust_vec = [800];
% number of seeds
num_seeds_vec = [200, 1000, 2000, 5000];
num_max_seeds = [200, 400, 600, 800]; %the maximum number of seeds to use in matching the clusters
num_params = length(num_seeds_vec);

% number of vertices
N = size(A,1);

% maximum dimension considered for the embedding
numdim = 30;

acc     = zeros(num_exp, num_params, num_runs);
runtime = zeros(num_exp, num_params, num_runs);
parfor r = 1:num_runs
  % temporary vectors for storing acc and runtime
  acc_		= zeros(num_exp, num_params);
  runtime_	= zeros(num_exp, num_params);
  
  for i = 1:num_params
	% make experiments reproducable
	rng(r*num_params+i);
	
	%load parameters
	% index vector for nonseed vertices
	num_seeds = num_seeds_vec(i);
	max_seeds = num_max_seeds(i);
	m = num_seeds;
	nonseeds = m+1:N;
	
	max_clust = max_clust_vec(1);
	
	[r max_clust]
	
	% random seeding
	Aperm = randperm(N);
	seedinds = Aperm(1:num_seeds);
	Bperm = randperm(N);
	% move vertices in B to match seeds in A
	for seed_ind = 1:num_seeds
		seed = seedinds(seed_ind);
		Bperm(Bperm==seed) = Bperm(seed_ind);
		Bperm(seed_ind) = seed;
	end
	Bend = Bperm(num_seeds+1:end);
	% create shuffled matrices
	AA = A(Aperm,Aperm);
	BB = B(Bperm,Bperm);
	
	% run experiments
	ex = 1;

	start = tic;
	% lsgm
	rng(r*num_params+i);

	[match clust_labels] = BigGMr_hardcoded( AA,BB,m, numdim, max_clust, max_seeds, 0,0,0);
	runtime_(ex,i) = toc(start);
	acc_(ex,i) = mean(Aperm(nonseeds)==Bperm(match(nonseeds)));
	ex = ex+1;
	
	% oracle accuracy
	Aperm_inv = 1:N;
	Aperm_inv(Aperm) = 1:N;
	Bperm_inv = 1:N;
	Bperm_inv(Bperm) = 1:N;
	acc_(ex,i) = mean(clust_labels(Aperm_inv(nonseeds),1)==clust_labels(Bperm_inv(nonseeds),2));
	ex = ex+1;
	[r max_clust mean(clust_labels(Aperm_inv(nonseeds),1)==clust_labels(Bperm_inv(nonseeds),2))]
    

  end
  runtime(:,:,r) = runtime_;
  acc(:,:,r) = acc_;
end

save(save_file_name, 'num_seeds_vec', 'num_max_seeds', 'numdim', 'max_clust_vec', 'acc', 'runtime');

