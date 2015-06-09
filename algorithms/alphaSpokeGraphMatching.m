function prob = alphaSpokeGraphMatching(A, B, s, topK, alpha, numrestarts, center, graphMatchAlg)
% prob = alphaSpokeGraphMatching(A, B, s, topK, alpha, numrestarts)
%
% Run seeded graph matching with alpha-spoke initialization to produce
% a one-to-many matching with estimated probabilities of match.
%
% INPUTS: A, B  : adjacency matrices of the two graphs to be matched
%         s     : number of seeds
%          topK : when TRUE return a matrix of potential matches, when
%                 FALSE return a vector of the best match
%         alpha : relative distance between center and permutation matrix
%   numrestarts : number of initialization restarts
%        center : the center of the alpha-spokes:
%                 'bari' uses the barycenter
%                 'convex' uses the solution to the convex relaxed graph
%                 matching problem
% graphMatchAlg : algorithm to use for graph matching
%
% OUTPUT: prob  : when topK is TRUE, a matrix whose i,j-th entry gives an 
%                 estimated probability that vertex i in graph A is matched
%                 to vertex j in graph B
%                 when topK is FALSE, a vector whose i-th entry is the 
%                 index of the best match in graph B for vertex i in 
%                 graph A

switch nargin
    case 4
        alpha = 0.2; % somewhat arbitrary default
        numrestarts = 100; % also arbitrary
        center = 'bari'; % default to barycenter
        graphMatchAlg = @seedgraphmatchell2; % SGM
    case 5
        numrestarts = 100; % also arbitrary
        center = 'bari'; % default to barycenter
        graphMatchAlg = @seedgraphmatchell2; % SGM
    case 6
        center = 'bari'; % default to barycenter
        graphMatchAlg = @seedgraphmatchell2; % SGM
    case 7
        graphMatchAlg = @seedgraphmatchell2; % SGM
end

[totv,~]=size(A); % number of vertices
n=totv-s; % number of non-seeds

id = eye(totv);
averageP = zeros(totv);

if (strcmp(center, 'bari'))
    centerMatrix = ones(n)/n;
elseif (strcmp(center, 'convex'))
    A22=A(s+1:s+n,s+1:s+n);
    B22=B(s+1:s+n,s+1:s+n);
    % This function doesn't use seeds correctly, so give it no seeds
    [~,centerMatrix]=relaxed_normAPPB_FW_seeds(A22,B22,0);
else
    centerMatrix = ones(n)/n;
end


parfor i = 1:numrestarts
    init = alphaSpokeInit(n, alpha, centerMatrix);
    [corr, ~] = graphMatchAlg(A, B, s, topK, init);
    %disp(sprintf('Finished matching iteratio %d', i));
    P = id(corr,:);
    averageP = averageP + P;
end
prob = averageP/numrestarts;

if (topK == false)
    % Return the best match
    [~, prob] = max(prob, [], 2);
end

end
    
