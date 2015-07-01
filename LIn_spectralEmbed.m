function [XA, XB] = Lin_spectralEmbed(A, B, numdim,embedMethod)
% this embeds two graphs into a common Euclidean space using adjacency
% spectral embedding
% numdim is the dimension to embed into
% A and B are the two graphs
if strcmp(embedMethod,'adjacency')
    [XA, sval, ~] = svds(A,numdim);
    XA = XA*sqrt(sval);
    [XB, sval, ~] = svds(B,numdim);
    XB = XB*sqrt(sval);
elseif strcmp(embedMethod,'laplacian')
    [XA, sval, ~] = svds(A,numdim*2,-0.000001);
    XA = XA(:,1:numdim)*sqrt(sval(1:numdim));
    [XB, sval, ~] = svds(B,numdim*2,-0.000001);
    XB = XB(:,1:numdim)*sqrt(sval(1:numdim));
else
    disp('4 parameters are expected')
end