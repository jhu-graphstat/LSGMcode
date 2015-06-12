function L = regularizedLaplacian(A)
% Compute the regularized laplacian of an adjacency matrix A\
%
% Reference: K. Chaudhuri, F. Chung, and A. Tsiatas. Spectral clustering of 
% graphs with general degrees in the extended planted partition model. 
% Journal of Machine Learning Research, pages 1–23,2012

rowSums = sum(A); % Row Sums
nrows = length(rowSums);
% Use average degree as regularizer: 
% Breaking the computation into two parts to help reduce possibility of
% round off errors for large graphs
regularizer = sum(rowSums/nrows)/(nrows - 1);
D_sqrt = sparse(1:nrows, 1:nrows, 1./sqrt(rowSums + regularizer)); % Sparse diagonal matrix

L = D_sqrt * A * D_sqrt;

end