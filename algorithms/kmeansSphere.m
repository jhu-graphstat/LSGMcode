function [IDX, centroid, Dis] = kmeansSphere(X, numclust, maxclustsize)
% Clusters the rows of X by first projecting onto the unit sphere and then
% applying kmeans clustering
%
% [IDX, centroid, Dis] = kmeansSphere(X, numclust)
%
% INPUT:
%           X : data matrix, each row is a separate data point
%     numclust: desired number of clusteres
%
% OUTPUT:
%         IDX : size(X,1)  vector of integers giving the cluster labels for
%               each data point
%    centroid : numclust x size(X,2) matrix, where each row i is the
%               centroid for cluster i
%         Dis : size(X,1) x numclust matrix, where entry i,j gives the
%               distance from data point i to cluster centroid j

startt = tic;
n = size(X,1); % number of data points
Y = zeros(size(X));
for i = 1:n
    Y(i,:) = X(i,:)/norm(X(i,:));
end

fprintf( 'done spherical projection: %f\n', toc(startt) );

[IDX, centroid, Dis] = kmeansHybrid(Y, numclust, maxclustsize);

end