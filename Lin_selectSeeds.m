function seeds = Lin_selectSeeds(pieceA,pieceB,s,s_max,ratio,num_ht)
% select seeds for SGM: active seed selection
seedsA = pieceA(s+1:end,1:s);
seedsB = pieceB(s+1:end,1:s);
avgdegA = mean(seedsA,1);
avgdegB = mean(seedsB,1);

%% compute hashtag seeds
s_max_ht = round(s_max*ratio);
seedcans = 1:num_ht;
 
% compute first seed
[val, ind] = min(abs(avgdegA(seedcans)-.5)+abs(avgdegB(seedcans)-.5));
seeds_ht = zeros(1,s_max_ht);
seeds_ht(1) = ind;
% set of remaining seed candidates
seedcans = setdiff(seedcans, ind);
p = [seedsA(:,ind); seedsB(:,ind)]/2;
for seed = 2:s_max_ht
    temp = repmat(p, [1,length(seedcans)]);  
    p_ = [seedsA(:,seedcans); seedsB(:,seedcans)]/seed+temp;
    p_1 = 1-p_;
    temp = p_.*log(p_) +p_1.*log(p_1);
    temp(find(isnan(temp))) = 0;
    H = -sum(temp,1);
    [val, ind] = max(H);
    p = p*seed/(seed+1);
    seeds_ht(seed) = seedcans(ind);
    seedcans = seedcans([1:ind-1, ind+1:end]);
end

%% compute agent seeds
s_max_at = s_max-s_max_ht;
seedcans = num_ht+1:s;

% compute first seed
[val, ind] = min(abs(avgdegA(seedcans)-.5)+abs(avgdegB(seedcans)-.5));
seeds_at = zeros(1,s_max_at);
seeds_at(1) = ind+num_ht;
% set of remaining seed candidates
seedcans = setdiff(seedcans, ind);
p = [seedsA(:,ind); seedsB(:,ind)]/2;
for seed = 2:s_max_at
    temp = repmat(p, [1,length(seedcans)]);  
    p_ = [seedsA(:,seedcans); seedsB(:,seedcans)]/seed+temp;
    p_1 = 1-p_;
    temp = p_.*log(p_) +p_1.*log(p_1);
    temp(find(isnan(temp))) = 0;
    H = -sum(temp,1);
    [val, ind] = max(H);
    p = p*seed/(seed+1);
    seeds_at(seed) = seedcans(ind);
    seedcans = seedcans([1:ind-1, ind+1:end]);
end
seeds = [seeds_ht,seeds_at];