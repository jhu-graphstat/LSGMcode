function seeds = activeSeedSelection(seedsA, seedsB, nSeeds)
% Based on the entropy of the connections in the adjacency submatrices
% seedsA and seedsB between seeds and nonseeds, activeSeedSelection finds a
% good collection of seeds.
%
% INPUTS:   seedsA : an s x m adjacency submatrix, where s = number of
%                    seeds and m = number of nonseeds
%           seedsB : an s x m adjacency submatrix, where s = number of
%                    seeds and m = number of nonseeds
%           nSeeds : number of seeds desired (nSeeds < s)
% 
% OUTPUTS:   seeds : a vector with indices giving the "best" seeds

s = size(seedsA,1);
seeds = NaN;

if size(seedsB,1) ~= s
    warning('Adjacency submatrices do not have the same number of seeds')
    return
end
% Possible improvements include: throwing away seeds with no connections
% instead of just halting when nothing is connected.
if max(seedsA(:)) == 0 && max(seedsB(:)) == 0
    warning('Nonseeds not connected to seeds, expect poor results: not using any seeds')
    return
end
avgdegA = mean(seedsA,2);
avgdegB = mean(seedsB,2);
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
p = [seedsA(ind,:), seedsB(ind,:)]/2;
for seed = 2:nSeeds
    temp = repmat(p, [length(seedcans),1]);
    p_ = [seedsA(seedcans,:), seedsB(seedcans,:)]/seed+temp;
    p_1 = 1-p_;
    temp = p_.*log(p_) +p_1.*log(p_1);
    temp(find(isnan(temp))) = 0;
    H = -sum(temp,2);

    [val, ind] = max(H);
    p = p*seed/(seed+1);

%			[val, ind] = min(abs(avgdegA(seedcans) -.5) +...
%							 abs(avgdegB(seedcans) -.5));

%			[val, ind] = min(abs(avgdegA(seedcans) -.5 +penaltyA) +...
%							 abs(avgdegB(seedcans) -.5 +penaltyB));
    seeds(seed) = seedcans(ind);
    seedcans = seedcans([1:ind-1, ind+1:end]);
end