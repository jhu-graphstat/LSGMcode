%% A function for Dan S to test that things kind of work

% Start in the main git repo directory
addpath(genpath('.'))

%% Load the data

% YOU NEED TO CHANGE THIS IF YOU ARE NOT DAN SUSSMAN

dir_04 = '~/Dropbox/Data/twitter_xdata_2015-05-16/2014_04/';
node_fn_04 = [dir_04 'twitter_prune_w50000.nodes.txt'];
edge_fn_04 = [dir_04 'twitter_prune_w50000.edges.txt'];
[node_04, adj_04] = readTwitter(node_fn_04,edge_fn_04);

dir_05 = '~/Dropbox/Data/twitter_xdata_2015-05-16/2014_05/';
node_fn_05 = [dir_05 'twitter_prune_w50000.nodes.txt'];
edge_fn_05 = [dir_05 'twitter_prune_w50000.edges.txt'];
[node_05, adj_05] = readTwitter(node_fn_05,edge_fn_05);

%% Get common vertices
[node_common, adj_04_common, adj_05_common] = ...
    commonVerticesTwitter(node_04,adj_04, node_05, adj_05);

%% Reduce to only users
node_common.is_hash = cellfun(@(name) name(1)=='#', node_common.name);

adj_04_user = adj_04_common(~node_common.is_hash,~node_common.is_hash);
adj_05_user = adj_05_common(~node_common.is_hash,~node_common.is_hash);
