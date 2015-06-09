function [ acc, runtime] = compareLSGM(num_runs, N, m, corrln, lam)
% Computes the accuracy and runtime of LSGM, SGM, U, RANK, QCV, rand, PATH
%Input: numsim: number of simulations
%   N: number of vertices
%   m: number of seeds
%   corrln: correlation
%   lam: block probability matrix
%Output: acc accuracy of 6 methods
%   runtime: computes the run time of 6 methods

%numdim = rank(lam);
numdim = 10;
max_clust_size = sum(N)/length(lam);
acc = zeros(2, num_runs);
runtime = zeros(2, num_runs);
for r = 1:num_runs
	% make experiments reproducable
	rng(r);
   
	% generate correlated graphs
	[A, B, shuffle] = sampleGraphs(m, N, corrln, lam);
	ex = 1;
	
	% lsgm with alpha-spoke top-K
	start = tic;
    alpha = 0.2;
    nRestarts = 10;
    topKMatching = @(A, B, m, topK) alphaSpokeGraphMatching(A, B, m, topK, ...
                               alpha, nRestarts, 'bari', @seedgraphmatchell2);
	matchMatrix = BigGMr( A,B,m,numdim, max_clust_size, @spectralEmbed, ...
                         @kmeansAlg, topKMatching, true);
    [weight, match] = max(matchMatrix, [], 2);
	runtime(ex,r)  = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)==match((m+1):end)');
	ex = ex+1;
    
    % lsgm w/o top-k
    start = tic;
   
	match = BigGMr( A,B,m,numdim, max_clust_size, @spectralEmbed, ...
                   @kmeansAlg, @seedgraphmatchell2, false);
  	runtime(ex,r)  = toc(start);
	acc(ex,r) = mean(shuffle(m+1:end)==match(m+1:end));
	ex = ex+1;
	
%	start = tic;
%	match = BigGM( A,B,m,N, numdim, max_clust_size, @spectralEmbed, @kmeansAlg, @seedgraphmatchell2test);
%	runtime(ex,r)  = toc(start);
%	acc(ex,r) = mean(shuffle(m+1:end)-m==match);
%	ex = ex+1;

end



end

