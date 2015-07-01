function [gA_,gB_] = embeddingKmean(A,B,s,max_clust_size,numdim)
indA = 1:length(A);
indB = 1:length(B);

% embed and cluster
[gA_,gB_] = EBKmean(A,B,indA,indB,s,max_clust_size,numdim); 
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

function [gA_,gB_] = EBKmean(A,B,indA,indB,s,max_clust_size,numdim) 
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
numclust = eva.OptimalK;
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
        num_seeds = length(seeds);
        gAaug = [seeds,gA_{i}];
        gBaug = [seeds,gB_{i}];
        AA = A(gAaug,gAaug);
        BB = B(gBaug,gBaug);
        [gA_{i},gB_{i}] = EBKmean(AA,BB,indA(gAaug),indB(gBaug),num_seeds,max_clust_size,numdim); 
    else
        gA_{i} = indA(gA_{i});
        gB_{i} = indB(gB_{i});
    end
end

