show_output = 1;


startt = tic;
nodeFn1 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_1_nodes.txt';
edgeFn1 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_1_edges_relabeled.txt';
nodeFn2 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_2_nodes.txt';
edgeFn2 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_2_edges_relabeled.txt';
twitterEdgeFn = '/home/jason/Documents/XDATA2015/summer-school/data/toLSGM_006_twitter_edges_october_relabeled.txt';

adj2= readTwitter2(edgeFn1);
adj1 = readTwitter2(twitterEdgeFn);
if show_output
	fprintf( 'done data injest: %f\n', toc(startt) );
end

startt = tic;
A = regularizedLaplacian(adj1);
B = regularizedLaplacian(adj2);
if show_output
	fprintf( 'done laplacian regularization: %f\n', toc(startt) );
end

numdim = 10;
s = 69;

nANonseeds = length(A) - s;
nBNonseeds = length(B) - s;

maxclustsize = 1000;
numclust = ceil(max(length(A), length(B))/maxclustsize);

startt = tic;
[XA XB] = spectralEmbed(A, B, numdim);
if show_output
	fprintf( 'done projection: %f\n', toc(startt) );
end

%% compute procrusties othogonal projection (on the seed vertices)
startt = tic;
[~,~,TRANSFORM]=procrustes(XA(1:s,:),XB(1:s,:));
TRANSFORM.c=ones(nBNonseeds+s,1)*TRANSFORM.c(1,:);
XB = TRANSFORM.b * XB * TRANSFORM.T + TRANSFORM.c;
if show_output
	fprintf( 'done procrusties: %f\n', toc(startt) );
end

%% cluster using the embedding
startt = tic;
XAXB=[XA;XB];
nonseedsA = s+1:s+nANonseeds;
nonseedsB = s+nANonseeds+ s+1:2*s + nANonseeds + nBNonseeds;
nonseeds = [nonseedsA, nonseedsB];
%[IDX, centroid, Dis] = clustAlg(XAXB, numclust);
[IDX, ~, Dis] = kmeansSphere(XAXB, numclust, maxclustsize);
if show_output
	fprintf( 'done clustering: %f\n', toc(startt) );
end


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
