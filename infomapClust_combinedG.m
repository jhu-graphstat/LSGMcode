function [gA_,gB_] = infomapClust_combinedG(A,B,s,ind_start,max_clust_size,trials,lambda)

num_seeds = s-ind_start; % number of non-hashtag seeds
AA = A(ind_start+1:end,ind_start+1:end);
BB = B(ind_start+1:end,ind_start+1:end);
indA = ind_start+1:length(A);
indB = indA;
% combined the two graphs AA and BB
A11 = AA(1:num_seeds,1:num_seeds);
A12 = AA(1:num_seeds,num_seeds+1:end);
A21 = AA(num_seeds+1:end,1:num_seeds);
A22 = AA(num_seeds+1:end, num_seeds+1:end);
B12 = BB(1:num_seeds,num_seeds+1:end);
B21 = BB(num_seeds+1:end,1:num_seeds);
B22 = BB(num_seeds+1:end, num_seeds+1:end);
Z = zeros(length(A22),length(B22));
G = [A11,A12,B12;A21,A22,Z;B21,Z',B22];

% match communities in the two graphs
ratio = 0.05;
[gA_,gB_] = matchCommunities(G,indA,indB,num_seeds,max_clust_size,ratio,trials,lambda);
% while sum(cellfun(@iscell,gA_))>0
%     gA_ = [gA_{:}];
%     gB_ = [gB_{:}];
% end
% 
% numclust = length(gA_);
% pieceA_ = cell(1,numclust);
% pieceB_ = cell(1,numclust);
% for i = 1:numclust
%     gAaug = [1:s, gA_{i}];
%     gBaug = [1:s, gB_{i}];
%     pieceA_{i}=A(gAaug,gAaug);
%     pieceB_{i}=B(gBaug,gBaug);
% end

function [dim,ind] = elbow(cellarray)
num = cellfun(@length,cellarray);
[tmp,ind] = sort(num,'descend');
% find the "elbow" in the curve
d1 = diff(tmp);
d2 = diff(d1);
[~, dim] = max(d2);
% find second elbow
d2(1:dim) = 0;
[~, dim2] = max(d2(dim:end));
dim = dim2+dim;
% find third elbow
d2(1:dim) = 0;
[~, dim2] = max(d2(dim:end));
dim = dim2+dim;

function cellarray = merge(cellarray,dim)
for i = dim+1:length(cellarray)
    cellarray{dim} = [cellarray{dim},cellarray{i}];
end
cellarray(dim+1:length(cellarray))=[];

function [gA_,gB_] = matchCommunities(G,indA,indB,num_seeds,max_clust_size,ratio,trials,lambda)
% run infomap on the combined graph
addpath('CommunityDetection')
com_tmp = infomap(G,trials,lambda);
numclust = max(com_tmp);
com = {};
com_seeds = com_tmp(1:num_seeds)';
com{1} = com_tmp(num_seeds+1:length(indA))';
com{2} = com_tmp(length(indA)+1:end)';
% assign vertices to communities
gA_ = cell(1,numclust);
gB_ = cell(1,numclust);
seeds_ = cell(1,numclust);
for i = 1:numclust
    gA_{i} = indA(find(com{1}==i)+num_seeds);
    gB_{i} = indB(find(com{2}==i)+num_seeds);
    seeds_{i} = find(com_seeds==i);
end
% find the large community clusters
[dim,ind] = elbow(seeds_);
seeds_ = merge(seeds_(ind),dim);
gA_ = merge(gA_(ind),dim);
gB_ = merge(gB_(ind),dim);

a = cellfun(@length,gA_)
b = cellfun(@length,gB_)

% run informap on the seeds
addpath('CommunityDetection')
Com_tmp = infomap(G(1:num_seeds,1:num_seeds),trials,1);
Com = cell(1,max(Com_tmp));
for i = 1:length(Com)
    Com{i} = find(Com_tmp==i)';
end
[dim,ind] = elbow(Com);
Com = Com(ind);
Com = merge(Com,dim);

% match/merge communities using seeds
n = length(Com)-1;
m = length(seeds_);
count = zeros(n,m);
for i = 1:n
    for j = 1:m
        count(i,j) = length(intersect(Com{i},seeds_{j}));
    end
end
count
count = diag(1./sum(count,2))*count


% % split up large communities
% numclust = length(gA_);
% tmp = max([cellfun(@length,gA_);cellfun(@length,gB_)]);
% pInd = find(tmp > max_clust_size);
% for i = 1:numclust
%     if sum(i == pInd)
%         num_seeds = [length(sA_{i}),length(sB_{i})];
%         gAaug = [sA_{i},gA_{i}];
%         gBaug = [sB_{i},gB_{i}];
%         AA = A(gAaug,gAaug);
%         BB = B(gBaug,gBaug);
%         [gA_{i},gB_{i}] = matchCommunities(AA,BB,indA(gAaug),indB(gBaug),num_seeds,max_clust_size,ratio,trials,lambda);
%     else
%         gA_{i} = indA(gA_{i});
%         gB_{i} = indB(gB_{i});
%     end
% end