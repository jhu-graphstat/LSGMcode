function seeds = selectSeeds(pieceA,pieceB,s,s_max)
% select seeds for SGM: active seed selection
seedsA = pieceA(s+1:end,1:s);
seedsB = pieceB(s+1:end,1:s);
avgdegA = mean(seedsA,1);
avgdegB = mean(seedsB,1);
seedcans = 1:s;
    
% compute first seed
[val, ind] = min(abs(avgdegA(seedcans)-.5)+abs(avgdegB(seedcans)-.5));
seeds = zeros(1,s_max);
seeds(1) = ind;
% set of remaining seed candidates
seedcans = setdiff(seedcans, ind);
p = [seedsA(:,ind); seedsB(:,ind)]/2;
for seed = 2:s_max
    temp = repmat(p, [1,length(seedcans)]);  
    p_ = [seedsA(:,seedcans); seedsB(:,seedcans)]/seed+temp;
    p_1 = 1-p_;
    temp = p_.*log(p_) +p_1.*log(p_1);
    temp(find(isnan(temp))) = 0;
    H = -sum(temp,1);
    [val, ind] = max(H);
    p = p*seed/(seed+1);
    seeds(seed) = seedcans(ind);
    seedcans = seedcans([1:ind-1, ind+1:end]);
end