function prob = alphaSpokeGraphMatching(A, B, s, alpha, numrestarts)
% function prob = alphaSpokeGraphMatching(A, B, s, alpha, numrestarts)
%
% Run seeded graph matching with alpha-spoke initialization to produce
% a one-to-many matching with estimated probabilities of match.
%
% INPUTS: A, B  : adjacency matrices of the two graphs to be matched
%         s     : number of seeds
%         alpha : relative distance between center and permutation matrix
%   numrestarts : number of initialization restarts
%
% OUTPUT: prob  : a matrix whose i,j-th entry gives an estimated
%                 probability that vertex i in graph A is matched to vertex
%                 j in graph B


[totv,~]=size(A); % number of vertices
n=totv-s; % number of non-seeds

id = eye(totv);
averageP = zeros(totv);

for i = 1:numrestarts
    init = alphaSpokeInit(n, alpha, 'bari');
    [corr, ~] = seedgraphmatchell2(A, B, s, init);
    P = id(corr,:);
    averageP = averageP + P;
end
prob = averageP/numrestarts;
