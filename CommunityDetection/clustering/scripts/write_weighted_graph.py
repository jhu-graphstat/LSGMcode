#!/usr/bin/env python

import argparse
import networkx as nx
import codecs

# Main driver: command line interface
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get basic subgraphs")
    parser.add_argument("--in_graph", help="input graph, gpickle format", type=str, required=True)
    parser.add_argument("--out_graph", help="output graph, Pajek format", type=str, required=True)
    parser.add_argument("--out_label", help="labels to deal with utf-8", type=str, required=True)
    args = parser.parse_args()
    in_graph_fn = args.in_graph
    out_graph_fn = args.out_graph
    out_label_fn = args.out_label

    g = nx.read_gpickle(in_graph_fn)

    # Write out label file
    out_label_file = codecs.open(out_label_fn, 'w', encoding='utf-8')
    node_num = 1
    for n in g.nodes_iter():
        g.node[n]['lbl'] = node_num
        out_label_file.write(u'{} {}\n'.format(node_num, n))
        node_num += 1
    out_label_file.close()

    # nx.write_pajek(g, out_graph_fn)

    out_graph_file = open(out_graph_fn, 'w')
    out_graph_file.write('*vertices {}\n'.format(node_num-1))
    for n in g.nodes_iter():
        out_graph_file.write('{} {}\n'.format(g.node[n]['lbl'], 'node-'+str(g.node[n]['lbl'])))
    out_graph_file.write('*arcs {}\n'.format(g.number_of_edges()))
    for e0, e1 in g.edges_iter():
        tot = 0.0
        for ky in g[e0][e1].keys():
            tot += g[e0][e1][ky]
        out_graph_file.write('{} {} {}\n'.format(g.node[e0]['lbl'], g.node[e1]['lbl'], tot))
    out_graph_file.close()
