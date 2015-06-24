function [node, adj] = readTwitter(nodeFn, edgeFn)
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
n_node = height(node)
node.idx = (1:n_node)';
node.id = node.id;
% A hash to access the node index from its id
hash = sparse(1,node.id+1,node.idx);
%% Load edges
edge = dlmread(edgeFn,' ');
edge(isnan(edge))= 0;
edge(:,1:2) = edge(:,1:2)+1;
edge(:,3) = sum(edge(:,3:5)')';

%% make adj

edgeOutIdx = full(hash(edge(:,1)))';
edgeInIdx = full(hash(edge(:,2)))';

adj = sparse(edgeOutIdx,edgeInIdx, edge(:,3));
end