function seeds = activeSeedSelection(A, B, nSeeds)
% Based on the entropy of the connections in the adjacency submatrices
% seedsA and seedsB between seeds and nonseeds, activeSeedSelection finds a
% good collection of seeds.
%
% INPUTS:   A : an s x m adjacency submatrix, where s = number of
%                    seeds and m = number of nonseeds
%           B : an s x m adjacency submatrix, where s = number of
%                    seeds and m = number of nonseeds
%           nSeeds : number of seeds desired (nSeeds < s)
% 
% OUTPUTS:   seeds : a vector with indices giving the "best" seeds

[s,nsA] = size(A);
[sB,nsB] = size(B);

if sB ~= s
    error('Adjacency submatrices do not have the same number of seeds')
end
% Possible improvements include: throwing away seeds with no connections
% instead of just halting when nothing is connected.
if max(A(:)) == 0 && max(B(:)) == 0
    warning('Nonseeds not connected to seeds, expect poor results: not using any seeds')
    seeds = [];
    return
end

avgdegA = mean(A,2);
avgdegB = mean(B,2);
seedcans = 1:s;
% compute first seed
[~, ind] = min(abs(avgdegA(seedcans) -.5) +...
                 abs(avgdegB(seedcans) -.5));

seeds = zeros(1,nSeeds);
seeds(1) = ind;

% set of remaining seed candidates
seedcans = setdiff(seedcans, ind);
%    	pA = seedsA(ind,:)/2;
%    	pB = seedsB(ind,:)/2;
p = [A(ind,:), B(ind,:)]/2;
for seed = 2:nSeeds
    temp = repmat(p, [length(seedcans),1]);
    p_ = [A(seedcans,:), B(seedcans,:)]/seed+temp;
    p_1 = 1-p_;
    temp = p_.*log(p_) +p_1.*log(p_1);
    temp(isnan(temp)) = 0;
    H = -sum(temp,2);

    [~, ind] = max(H);
    p = p*seed/(seed+1);

%			[val, ind] = min(abs(avgdegA(seedcans) -.5) +...
%							 abs(avgdegB(seedcans) -.5));

%			[val, ind] = min(abs(avgdegA(seedcans) -.5 +penaltyA) +...
%							 abs(avgdegB(seedcans) -.5 +penaltyB));
    seeds(seed) = seedcans(ind);
    seedcans = seedcans([1:ind-1, ind+1:end]);
end
%		seeds = ind(1:num_seeds);
        