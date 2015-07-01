function [XA, XB] = spectralEmbed(A, B, numdim)
% this embeds two graphs into a common Euclidean space using adjacency
% spectral embedding
% numdim is the dimension to embed into
% A and B are the two graphs
[XA, sval, ~] = svds(A,numdim);
XA = XA*sqrt(sval);
[XB, sval, ~] = svds(B,numdim);
XB = XB*sqrt(sval);