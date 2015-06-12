function [IDX, centroid, Dis] = kmeansHybrid(X, numclust, maxclustsize)

[IDX, centroid] = kmeansBisect(X', numclust, maxclustsize);
IDX = IDX';
centroid = centroid';
Dis = zeros(size(X,1), size(centroid,1));
for c = 1:numclust
	Dis(:,c) = sum( (X -repmat(centroid(c,:), [size(X,1),1])).^2, 2 );
end
