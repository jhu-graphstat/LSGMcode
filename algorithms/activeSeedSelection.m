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

A = uint64(A>0);
B = uint64(A>0);
seeds = zeros(1,nSeeds);

% First find the the seed that is closest to being connected to half of
% the non-seeds in both graphs, ie highest entropy seed
avgdegA = mean(A,2);
avgdegB = mean(B,2);
seedCand = 1:s;
% compute first seed
[~, nextSeed] = min(abs(avgdegA(seedCand) -.5) +...
                 abs(avgdegB(seedCand) -.5));
seeds(1) = nextSeed;

% Remove from candidate seeds
seedCand = setdiff(seedCand, nextSeed);

% If we are looking for lots of seeds we will only do smart things for a 
% little while 
if nSeeds > 64
    initSeeds = 64;
    doAfter = true;
else
    initSeeds = nSeeds;
    doAfter = false;
end
    

pA = A(nextSeed,:);
pB = B(nextSeed,:);

% each new seed will be the seed that increases the entropy of the matrices 
% A(seeds,:) and B(seeds,:) viewed as a collection of column vectors
for seed = 2:initSeeds
    pA = bitshift(pA,1);
    pB = bitshift(pB,1);
    % Convert the columns into integers
    tempA = repmat(pA,[length(seedCand),1]);
    tempA = bitor(A(seedCand,:),tempA);
    
    tempB = repmat(pB,[length(seedCand),1]);
    tempB = bitor(A(seedCand,:),tempB);
    
    % Compute the entropy of each new seed and pick the best one
    [~,nextSeed] = max(rowEntropy(tempA)+rowEntropy(tempB));
    seeds(seed) = seedCand(nextSeed);
    seedCand = setdiff(seedCand,seedCand(nextSeed));
    
    % Add on to pA and pB for the next step
    pA = bitor(A(nextSeed,:),pA);
    pB = bitor(B(nextSeed,:),pB);
    if numel(unique(pA))==nsA || numel(unique(pB))==nsB && seed < initSeeds
        disp(seed)
        doAfter = true;
        break
    end
end

if doAfter
    % If there are lots of seeds then for the rest we'll just try to
    % maximize the sum of the entropies of each column
    A = double(A);
    B = double(B);
    p = mean([A(seeds(1:seed),:),B(seeds(1:seed),:)]);
    seedStart = seed+1;
    for seed = seedStart:nSeeds
        temp = repmat(p, [length(seedCand),1]);
        p_ = [A(seedCand,:), B(seedCand,:)]/seed+temp;
        p_1 = 1-p_;
        temp = p_.*log(p_) +p_1.*log(p_1);
        temp(isnan(temp)) = 0;
        H = -sum(temp,2);

        [~, nextSeed] = max(H);
        p = p*seed/(seed+1);
        seeds(seed) = seedCand(nextSeed);
        seedCand = seedCand([1:nextSeed-1, nextSeed+1:end]);
    end
end

end

function re = rowEntropy(A)
    [nr,nc] = size(A);
    re = zeros(nr,1);
    for r=1:nr
        [~,~,iu] = unique(A(r,:));
        freq = zeros(max(iu),1);
        for c=1:nc
            freq(iu(c)) = freq(iu(c))+1;
        end
        freq = freq/nc;
        re(r) = -sum(freq.* log(freq));
    end
    
end