function A = reweight(A,weighting)
% Reweights the nonzero elements of the matrix A
% 
% Anew = reweight(A,weighting)
% 
% Anew is the new matrix reweighted.
% 
% A is the original sparse or dense matrix
% weighting is a function handle or the string 'rank'. The non-zero
% elements of A are transformed according to the funciton. If weighting is
% equal to 'rank' then this the function @(s)tiedrank(s)/numel(s)
if nargin < 2
    weighting = 'rank';
end

[i,j,s] = find(A);

if ~iscell(weighting)
    weighting = {weighting};
end

for weight=weighting
    w = weight{1};
    if isa(w, 'function_handle')
        s = w(s);
    elseif strcmp(w,'rank')
        s = tiedrank(s)*2; % times two to get integers 
    end
end

A(sub2ind(size(A),i,j)) = s;

end