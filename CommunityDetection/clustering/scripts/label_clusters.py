#!/usr/bin/env python

import networkx as nx
import re

tgt_fn = 'info/targets.key'
graph_in_fn = 'graphs/sn_eg_edge_networkx.txt'
cluster_in_fn = 'tmp/clusters.txt'
graph_out_fn = 'graphs/sn_eg.gexf'

# Read in clusters
clustfile = open(cluster_in_fn, 'r')
clust_map = {}
for ln in clustfile:
    ln = ln.rstrip()
    f = ln.split()
    clust_map[f[0]] = f[1]
clustfile.close()

# Read in graph and label it
G = nx.read_edgelist(graph_in_fn)
not_found = 0
found = 0
for node in G.nodes():
    # Put in cluster num
    if (clust_map.has_key(node)):
        G.node[node]['cluster_num'] = clust_map[node]
    else:
        print "this shouldn't happen, cluster not found: {0}".format(node)
        raise

# print 'found : {0} / {1}'.format(found, found+not_found)

# Write out result
nx.write_gexf(G, graph_out_fn)
