function init = alphaSpokeInit(k, alpha, center)
% function init = alphaSpokeInit(k, alpha, center)
%
% Returns the alpha-spoke random initialization:
% (1-alpha)P + (alpha)C
% where P is a (uniformly) random permutation matrix, and C is the central
% matrix, either the baricenter or the solution to the convex relaxation
%

id = eye(k);
P = id(randperm(k),:);
if (strcmp(center, 'bari'))
    C = ones(k)/k;
else
    C = ones(k)/k;
end

init = (1-alpha)*P + alpha*C;