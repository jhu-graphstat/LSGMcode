function [nodeCommon, adj1,adj2] = ...
    commonVerticesTwitter(node1, adj1, node2, adj2, key)

% COMMONVERTICESTWITTER Gets the graphs corresponding to the vertex
% interesection of of two twitter grahps
%   Use the COMMONVERTICESTWITTER function to get the induced subgraphs
%   corresponding to the common vertices between two twitter graphs.
%
%   [NODECOMMON, ADJ1,ADJ2] = COMMONVERTICESTWITTER(NODE1, ADJ1, NODE2, ADJ2)
%   returns the induced subgraphs corresponding to the common nodes of the
%   two graphs. It performs an inner join on node1 and node2 based on the
%   'name' columns of the two tables and then collects the subgraphs
%   corresponding to the join based on the 'idx' columns.


if nargin < 5
    key = 'name';
end

nodeCommon = innerjoin(node1,node2,'Keys',key,...
    'LeftVariables',{key,'idx'},'RightVariables',{'idx'});
nodeCommon.is_hash = cellfun(@(name) name(1)=='#', nodeCommon.name);

adj1= adj1(nodeCommon.idx_node_1,nodeCommon.idx_node_1);
adj2 = adj2(nodeCommon.idx_node_2,nodeCommon.idx_node_2);


end

