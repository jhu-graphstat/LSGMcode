function adj = readTwitter2(edgeFn)
% READTWITTER Read twitter data from two text files 
%   Use the READTWITTER function to read in twitter data from a text file
%   for the nodes and another for the edges. READTWITTER returns a table
%   for the nodes and a sparse matrix for the adjacency matrix.
%
%   ADJ = READTWITTER(EDGEFN) creates a table and a sparse
%   adjacency matrix corresponding to the data in the files specified.

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

% Make a symmetric adjacency matrix
adj = sparse([edge.out, edge.in],[edge.in, edge.out], 1);
end