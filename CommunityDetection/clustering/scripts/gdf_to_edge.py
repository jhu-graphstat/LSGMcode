#!/usr/bin/env python

# import networkx as nx
import re

graph_in_fn = 'graphs/sn_eg.gdf'
graph_out_fn = 'graphs/sn_eg_edge_networkx.txt'
graph_out2_fn = 'graphs/sn_eg_edge.txt'

graph_in = open(graph_in_fn, 'r')
graph_out = open(graph_out_fn, 'w')
graph_out2 = open(graph_out2_fn, 'w')
edges = False
for ln in graph_in:
    ln = ln.rstrip()
    if (edges):
        f = ln.split(',')
        # graph_out.write('{0}'.format(ln) + ',{\'weight\': 1.0}\n')
        graph_out.write('{0} {1}\n'.format(f[0],f[1]))
        graph_out2.write('{0}\n'.format(ln))
    if (re.match("^edgedef", ln)):
        edges = True
graph_in.close()
graph_out.close()

