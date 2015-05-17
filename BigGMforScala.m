function [ match, clust_labels ] = BigGMforScala( A, B, s, numdim, max_clust_size)


start = tic;
sumn = length(A)-s;
numclust = ceil(sumn/max_clust_size);

% maxmium number of seeds to use
s_max = min(200,s);

% show output
show_output = false;

% perform embedding
startt = tic;
[XA XB] = spectralEmbed(A, B, numdim);
%if show_output
	fprintf( 'done projection: %f\n', toc(startt) );
%end

% compute procrusties othogonal projection (on the seed vertices)
startt = tic;
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(sumn+s,1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
%if show_output
	fprintf( 'done procrusties: %f\n', toc(startt) );
%end


% cluster using the embedding
startt = tic;
XAXB=[XA;XB];
nonseedsA = s+1:s+sumn;
nonseedsB = s+sumn+ s+1:2*(s+sumn);
[IDX, centroid, Dis] = kmeans0(XAXB, numclust);
fprintf( 'done clustering: %f\n', toc(startt) );

%% fix cluster sizes to be equal in both graphs

[pieceA_,pieceB_,gA_,gB_] = fixClusterSize(IDX, Dis, numclust,nonseedsA, nonseedsB)
clear IDX Dis

%% perform graph matching in parallel
startt = tic;
match = zeros(s+sumn,numclust);
tic
parfor i = 1:numclust
%parfor i = 1:numclust
	% load subgraph adjacency matrix
	gA = gA_{i};
	gB = gB_{i};
	pieceA = pieceA_{i};
	pieceB = pieceB_{i};
	
	% if graph is empty
	if length(gA)==0
		continue;
	end
	
    % perform graph match
    %tic
	ord = seedgraphmatchell2(pieceA, pieceB, s);
	%toc
    ord = ord(s+1:end)-s;
	
	% save results
	temp = zeros(s+sumn,1);
	temp(gA) = gB(ord);
    match(:,i) = temp;
end
toc
% combine results
match(1:s,1) = 1:s;
match = sum(match,2)';

fprintf( 'done matching: %f\n', toc(startt) );

