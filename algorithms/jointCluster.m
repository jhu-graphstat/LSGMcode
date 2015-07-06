function [IDX, nonseedsA, nonseedsB, XAXB] = jointCluster(A, B, s, nDim, maxClustSize)


nANonseeds = length(A) - s;
nBNonseeds = length(B) - s;

sumn = max(nANonseeds, nBNonseeds);

% number of clusters determined by max cluster size
numclust = ceil(sumn/maxClustSize);
% maxmium number of seeds to use
s_max = min(200,s);
% show output
 show_output = true;
 
startt = tic;

[XA, XB] = spectralEmbed(A, B, nDim);
if show_output
	fprintf( 'done projection: %f\n', toc(startt) );
end

%% compute procrusties othogonal projection (on the seed vertices)
if s > 0
    startt = tic;
    [~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
    TRANSFORM.c=ones(nBNonseeds+s,1)*TRANSFORM.c(1,:);
    XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
    if show_output
        fprintf( 'done procrusties: %f\n', toc(startt) );
    end
end

%% cluster using the embedding
startt = tic;
XAXB=[XA;XB];
nonseedsA = s+1:s+nANonseeds;
nonseedsB = s+nANonseeds+ s+1:2*s + nANonseeds + nBNonseeds;
nonseeds = [nonseedsA, nonseedsB];
%[IDX, centroid, Dis] = clustAlg(XAXB, numclust);
%IDX = kmeansSphere(XAXB(nonseeds,:), numclust, maxClustSize
IDX = kmeansSphere(XAXB, numclust, maxClustSize);

if show_output
    fprintf('done initial clustering: %f\n', toc(startt))
end

startt = tic;

% If any clusters are too large, we recursively cluster them
c = 1;
% Find all indicices belonging to cluster c
clusterC = find(IDX == c);
% Only keep non seeds from graph A and B
clusterC = clusterC(clusterC > s & (clusterC <= (s + nANonseeds) | clusterC > (2*s + nANonseeds)));
while(~isempty(clusterC))
    iA = clusterC(clusterC <= (s + nANonseeds) & clusterC > s);
    iB = clusterC(clusterC > nANonseeds + 2*s) - (nANonseeds + 2*s);
    nA = length(iA);
    nB = length(iB);
    if  (nA > 0 && nB > 0) && (length(iA) + length(iB) > 2*maxClustSize)
        %numNewClusters = ceil(length(clusterC)/maxClustSize);
        seeds = activeSeedSelection(A(1:s, iA), B(1:s, iB), s_max);
        %seeds = 1:s; % temporarily use all seeds
        if ~isempty(seeds) 
            iA = [seeds'; iA];
            iB = [seeds'; iB];
        else
            s_max = 0; % No useful seeds
        end
        [subIDX, ~, ~] = jointCluster(A(iA,iA), B(iB,iB), s_max, nDim, maxClustSize);
        % Only grab nonseeds
        subIDX = subIDX([(s_max+1):(s_max+nA), (2*s_max+nA + 1):(2*s_max+nA+nB)]);
        numNewClusters = max(subIDX);
        % Shift old cluster labels
        IDX(IDX > c) = IDX(IDX > c) + numNewClusters - 1;
        % Assign labels c, c+1, ... to the newly split cluster
        IDX(clusterC) = subIDX + c - 1;
        c = c + numNewClusters - 1; % skip to a new cluster
    end
    c = c + 1;
    clusterC = find(IDX == c);
    clusterC = clusterC(clusterC > s & (clusterC <= (s + nANonseeds) | clusterC > (2*s + nANonseeds)));
end

%nonseedsA = 1:nANonseeds;
%nonseedsB = (nANonseeds + 1):(nANonseeds + nBNonseeds);

if show_output
    fprintf('done recursive clustering: %f\n', toc(startt))
end
end

