function idx = getLargestCC(A)
[ci, size] = components(A);
[~, index] = max(size);
idx = find(ci == index);
end