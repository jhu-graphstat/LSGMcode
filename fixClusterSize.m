function [pieceA_,pieceB_,gA_,gB_] = fixClusterSize(A,B,IDX, Dis, numclust,nonseedsA, nonseedsB)

sumn = length(nonseedsA);
s = length(A)-sumn;

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

clustsizes=round((clustsizesA+clustsizesB)/2);

% add or subtract 1 from largest cluster sizes
[val,ind] = sort(clustsizes, 1, 'descend');
temp = sumn-sum(clustsizes);
mask = ind(1:abs(temp));
clustsizes(mask) = clustsizes(mask)+sign(temp);

% save clusters(subgraphs) to match (so it can be parallelized)
pieceA_ = cell(1,numclust);
pieceB_ = cell(1,numclust);
gA_ = cell(1,numclust);
gB_ = cell(1,numclust);

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
	
    gAaug = [1:s, (gA+s)'];
    gBaug = [1:s, (gB+s)'];
    
    pieceA=A(gAaug,gAaug);
    pieceB=B(gBaug,gBaug);
    
    pieceA_{i} = pieceA;
    pieceB_{i} = pieceB;
    
    gA_{i} = gA;
    gB_{i} = gB;
end