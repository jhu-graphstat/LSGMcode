function [gA_,gB_] = infomapClust(A,B,s,max_clust_size,numdim,trials,lambda)
indA = 1:length(A);
indB = 1:length(B);
% match communities of the two graphs
ratio = 0.01;
[gA_,gB_] = matchCommunities(A,B,indA,indB,s,max_clust_size,ratio,numdim,trials,lambda);
% un-nest gA_ and gB_
while sum(cellfun(@iscell,gA_))>0
    ind = find(cellfun(@iscell,gA_));
    for i = ind
        tmp = {gA_,gA_{i}};
        gA_ = [tmp{:}];
        tmp = {gB_,gB_{i}};
        gB_ = [tmp{:}];
    end
    gA_(ind) = [];
    gB_(ind) = [];
end
cellfun(@length,gA_)
cellfun(@length,gB_)

gA_
gB_
gA = [];
gB = [];
for i = 1:length(gA_)
    gA = [gA,gA_{i}];
    gB = [gB,gB_{i}];
end
plot(sort(gA))
hold on
plot(sort(gB),'r')


function [seeds,non_seeds] = findInd(com,num_seeds)
numclust = max(com);
seeds = cell(1,numclust);
non_seeds = cell(1,numclust);
for i = 1:numclust
    tmp = find(com == i);
    seeds{i} = tmp(tmp <= num_seeds)';
    non_seeds{i} = tmp(tmp > num_seeds)';
end

function [dim,ind] = elbow(num)
[tmp,ind] = sort(num,'descend');
% find the "elbow" in the curve
d1 = diff(tmp);
d2 = diff(d1);
[~, dim] = max(d2);
for i = 1:3
    d2(1:dim) = 0;
    [~, dim2] = max(d2(dim:end));
    dim2 = dim2+dim-1;
    if tmp(dim2) > 1
        dim = dim2;
    else
        continue
    end
end


function [gA_,gB_] = matchCommunities(A,B,indA,indB,num_seeds,max_clust_size,ratio,numdim,trials,lambda)
% run infomap algorithm
com = {};
addpath('CommunityDetection')
com{1} = infomap(A,trials,lambda);
com{2} = infomap(B,trials,lambda);

% generate sets of seed indices
[seedsA,comA] = findInd(com{1},num_seeds);
[seedsB,comB] = findInd(com{2},num_seeds);

% match communities using seeds
n = length(seedsA);
m = length(seedsB);
count = zeros(n,m);
for i = 1:n
    for j = 1:m
        count(i,j) = length(intersect(seedsA{i},seedsB{j}));
    end
end

% find the large community clusters
% find the "elbow" in the curve
[tmp,indr] = sort(cellfun(@length,seedsA),'descend');
[tmp,indc] = sort(cellfun(@length,seedsB),'descend');

% match clusters
count = count(indr,indc);
if n>m
    count(m,:) = sum(count(m:end,:));
    count(m+1:end,:) = [];
elseif n<m
    count(:,n) = sum(count(:,n:end),2);
    count(:,n+1:end) = [];
end
corr_c = lapjv(-count', 0.01);
count = count(corr_c,:);
if length(count) <= 2
    disp('numclust = 1 (not enough seeds)')
    numclust = sum(diag(count)>0);
else
    [dim,ind] = elbow(diag(count));
    numclust = dim+1;
end
% update count matrix
count(numclust,:) = sum(count(numclust:end,:));
count(numclust+1:end,:) = [];
count(:,numclust) = sum(count(:,numclust:end),2);
count(:,numclust+1:end) = [];

% assign vertices to clusters
gA_ = cell(1,numclust);
gB_ = cell(1,numclust);
sA_ = cell(1,numclust);
sB_ = cell(1,numclust);
for i = 1:numclust-1
    gA_{i} = comA{indr(corr_c(i))};
    gB_{i} = comB{indc(i)};
    sA_{i} = seedsA{indr(corr_c(i))};
    sB_{i} = seedsB{indc(i)};
end
tmp = setdiff(1:max(com{1}),indr(corr_c(1:numclust-1)));
gA = [];
sA = [];
for i = 1:length(tmp)
    gA = [gA,comA{tmp(i)}];
    sA = [sA,seedsA{tmp(i)}];
end
gA_{numclust} = gA;
sA_{numclust} = sA;

tmp = setdiff(1:max(com{2}),indc(1:numclust-1));
gB = [];
sB = [];
for i = 1:length(tmp)
    gB = [gB,comB{tmp(i)}];
    sB = [sB,seedsB{tmp(i)}];
end
gB_{numclust} = gB;
sB_{numclust} = sB;

% merging communities 
flag = 1;
threshold = ceil(num_seeds*ratio);
while flag
    ctmp = count;
    d = diag(ctmp)';
    ctmp = ctmp - diag(d); % set diagonal to zero
    [val,c1] = max(ctmp);
    [val,c2] = max(val);
    c1 = c1(c2);
    if (ctmp(c1,c2) > threshold)&&(length(d) > 1)
        % update the count matrix 
        count(:,c1) = count(:,c1)+count(:,c2);
        count(:,c2) = [];
        count(c1,:) = count(c1,:)+count(c2,:);
        count(c2,:) = [];
        % merging the two clusters
        gA_{c1} = [gA_{c1},gA_{c2}];
        gA_(c2) = [];
        gB_{c1} = [gB_{c1},gB_{c2}];
        gB_(c2) = [];
        sA_{c1} = [sA_{c1},sA_{c2}];
        sA_(c2) = [];
        sB_{c1} = [sB_{c1},sB_{c2}];
        sB_(c2) = [];
    else
        flag = 0;
    end
end

% split up large communities
numclust = length(gA_);
tmp = max([cellfun(@length,gA_);cellfun(@length,gB_)]);
pInd = find(tmp > max_clust_size);
for i = 1:numclust % Due to memory issue, for loop is used instead of parfor for large networks
    seeds = union(sA_{i},sB_{i});
    if sum(i == pInd) && length(seeds)>3
        if numclust == 1 % when the community detection algorithm can no longer break up the networks
            % use embed and kmean 
            disp('...kmean algorithm...')
            num_seeds = length(seeds);
            gAaug = [seeds,gA_{i}];
            gBaug = [seeds,gB_{i}];
            AA = A(gAaug,gAaug);
            BB = B(gBaug,gBaug);
            [gA_{i},gB_{i}] = EBKmean(AA,BB,indA(gAaug),indB(gBaug),num_seeds,max_clust_size,ratio,numdim,trials,lambda);
        else
            % infomapClust
            num_seeds = length(seeds);
            gAaug = [seeds,gA_{i}];
            gBaug = [seeds,gB_{i}];
            AA = A(gAaug,gAaug);
            BB = B(gBaug,gBaug);
            ratio = 1.2*ratio;
            [gA_{i},gB_{i}] = matchCommunities(AA,BB,indA(gAaug),indB(gBaug),num_seeds,max_clust_size,ratio,numdim,trials,lambda);
        end   
    else
        gA_{i} = indA(gA_{i});
        gB_{i} = indB(gB_{i});
    end
end

function [gA_,gB_] = EBKmean(A,B,indA,indB,s,max_clust_size,ratio,numdim,trials,lambda) 
% embedding
[XA, XB] = spectralEmbedElbow_directed(A, B,numdim,0);

% compute procrusties othogonal projection (on the seed vertices)
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(length(B),1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;

% cluster (K-means) using the embedding
XAXB=[XA;XB];
nonseedsA = s+1:length(A);
nonseedsB = (length(A)+s+1):(length(A)+length(B));
nonseeds = [nonseedsA, nonseedsB];
eva = evalclusters(XAXB(nonseeds,:),'kmeans','DaviesBouldin','KList',[1:10]);
numclust = eva.OptimalK
[IDX, centroid, Dis] = kmeansAlgr(XAXB, numclust, nonseeds);
IDXA = IDX(s+1:length(A))';
IDXB = IDX(length(A)+s+1:end)';
IDSA = IDX(1:s)';
IDSB = IDX(length(A)+1:length(A)+s)';
gA_ = cell(1,numclust);
gB_ = cell(1,numclust);
sA_ = cell(1,numclust);
sB_ = cell(1,numclust);
for i=1:numclust
   gA_{i} = find(IDXA==i)+s;
   gB_{i} = find(IDXB==i)+s;
   sA_{i} = find(IDSA==i);
   sB_{i} = find(IDSB==i);
end

% split up large communities
numclust = length(gA_);
tmp = max([cellfun(@length,gA_);cellfun(@length,gB_)]);
pInd = find(tmp > max_clust_size);
for i = 1:numclust % Due to memory issue, for loop is used instead of parfor for large networks
    seeds = union(sA_{i},sB_{i});
    if sum(i == pInd) && length(seeds)>3
        % infomapClust
        num_seeds = length(seeds);
        gAaug = [seeds,gA_{i}];
        gBaug = [seeds,gB_{i}];
        AA = A(gAaug,gAaug);
        BB = B(gBaug,gBaug);
        ratio = 1.5*ratio;
        [gA_{i},gB_{i}] = matchCommunities(AA,BB,indA(gAaug),indB(gBaug),num_seeds,max_clust_size,ratio,numdim,trials,lambda);
    else
        gA_{i} = indA(gA_{i});
        gB_{i} = indB(gB_{i});
    end
end
