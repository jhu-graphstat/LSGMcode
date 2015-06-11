function [node, adj] = readTwitter2(nodeFn, edgeFn)
% READTWITTER Read twitter data from two text files 
%   Use the READTWITTER function to read in twitter data from a text file
%   for the nodes and another for the edges. READTWITTER returns a table
%   for the nodes and a sparse matrix for the adjacency matrix.
%
%   [NODE, ADJ] = READTWITTER(NODEFN, EDGEFN) creates a table and a sparse
%   adjacency matrix corresponding to the data in the files specified.

%% Load nodes
node = readtable(nodeFn,'ReadVariableNames',false,'Delimiter',' ');
node.Properties.VariableNames = {'id' 'name'};
n_node = height(node);
node.idx = (1:n_node)';
% A hash to access the node index from its id
% Add 1 since Kitware data is initialized at zero
hash = containers.Map(node.id, node.idx);
%% Load edges
edge = readtable(edgeFn,'ReadVariableNames',false,'Delimiter',' ');

%% Make nice
edge = edge(:,1:2); % Kitware data doesn't have attributes
edge.Properties.VariableNames = {'out','in'};

%% make adj

%edge.out_idx = full(hash(edge.out))'; % remember node list init at zero
%edge.in_idx = full(hash(edge.in))';

if length(edge.out) ~= length(edge.in)
    error('Number of out vertices does not match in vertices')
end

nEdges = length(edge.out);

edge.out_idx = zeros(nEdges,1);
edge.in_idx = zeros(nEdges,1);


for i = 1:nEdges
    edge.out_idx(i) = hash(edge.out(i));
    edge.in_idx(i) = hash(edge.in(i));
end

% Make a symmetric adjacency matrix
adj = sparse([edge.out_idx, edge.in_idx],[edge.in_idx, edge.out_idx], 1);
end