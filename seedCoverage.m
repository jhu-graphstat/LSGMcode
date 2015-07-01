function [pA,pB] = seedCoverage(seedinds,A,B)
num_seeds = length(seedinds);
aC = zeros(1,length(A));
aC(seedinds) = 1;
bC = zeros(1,length(B));
bC(seedinds) = 1;
for i = 1:num_seeds
    ind = seedinds(i);
    % find one-hop neighbors 
    aC(find(A(ind,:)))=1;
    bC(find(B(ind,:)))=1;
end
pA = sum(aC)/size(A,1);
pB = sum(bC)/size(B,1);