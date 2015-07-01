function [XA, XB] = spectralEmbedElbow_directed(A, B, numdim, projectToSphere)

if nargin < 4
    warning('projectToSphere not specified, using true as default');
    projectToSphere = true;
end

%% Elbows for A

% embed A to numDim+2 dimension
[XA_l, svalA_l, ~] = svds(A, min([ceil(numdim/2)+1, length(A)]) );
[XA_r, svalA_r, ~] = svds(A', min([ceil(numdim/2)+1, length(A)]) );
svalA = 0.5*(svalA_l+svalA_r);

% find first and second derivatives
d1_A = diff(diag(svalA));
d2_A = diff(d1_A)';

% find the first "elbow" in the curve via max of second differences
[~, dimA]  = max(d2_A);
% find the second elbow as the next highest after the first
d2_A(1:dimA) = 0;
[~, dimA] = max(d2_A);

%% Elbows for B

% embed B
[XB_l, svalB_l, ~] = svds(B, min([ceil(numdim/2)+1, length(B)]) );
[XB_r, svalB_r, ~] = svds(B', min([ceil(numdim/2)+1, length(B)]) );
svalB = 0.5*(svalB_l+svalB_r);

% find the "elbow" in the curve
d1_B = diff(diag(svalB));
d2_B = diff(d1_B)';
[~, dimB] = max(d2_B);
% find second elbow
d2_B(1:dimB) = 0;
[~, dimB] = max(d2_B(dimB:end));

%% use max of 2 elbows
dim = max(dimA, dimB);


%% Scale the eigenvectors to get the embeddings for both graphs
XA = [XA_l(:,1:dim)*sqrt(svalA(1:dim,1:dim)),XA_r(:,1:dim)*sqrt(svalA(1:dim,1:dim))];
XB = [XB_l(:,1:dim)*sqrt(svalB(1:dim,1:dim)),XB_r(:,1:dim)*sqrt(svalB(1:dim,1:dim))];

if(projectToSphere)
    XA = XA *diag(sum(XA.^2).^(-.5));
    XB = XB *diag(sum(XB.^2).^(-.5));
end