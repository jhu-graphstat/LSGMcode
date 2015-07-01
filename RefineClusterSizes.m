function [nonseedsA,nonseedsB,flag] = RefineClusterSizes(IDX, Dis, numclust,nonseedsA, nonseedsB,s)
flag = 0;
% fixing cluster sizes to be equal in both graphs
IDXA = IDX(nonseedsA);
DisA = Dis(nonseedsA,:);
IDXB = IDX(nonseedsB);
DisB = Dis(nonseedsB,:);
clear IDX Dis

clustsizesA=zeros(numclust,1);
clustsizesB=zeros(numclust,1);
for i=1:numclust
   iiAi=find(IDXA==i);
   iiBi=find(IDXB==i);
   clustsizesA(i)=length(iiAi);
   clustsizesB(i)=length(iiBi);
end

if mean(clustsizesA==clust) == 1
    flag = 1;
else
    % set the cluster sizes
    clustsizes = min(clustsizesA,clustsizesB);

    % save nodes to match 
    nodesA = [];
    nodesB = [];

    for i=1:numclust
        if clustsizes(i) == 0
            continue
        end
        % sort by distance
        [~, TA] = sort(DisA(:,i));
        [~, TB] = sort(DisB(:,i));
        % use closest vertices
        gA = TA(1:clustsizes(i));
        gB = TB(1:clustsizes(i));
        % make these vertices unselectable for next time
        DisA(gA,:) = inf;
        DisB(gB,:) = inf;
        nodesA = [nodesA,gA+s];
        nodesB = [nodesB,gB+s];
    end
    nonseedsA = sort(nodesA);
    nonseedsB = sort(nodesB);
end