#!/usr/bin/env python

import argparse
import networkx as nx
import codecs

# Main driver: command line interface
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get basic subgraphs")
    parser.add_argument("--in_map", help="input module map", type=str, required=True)
    parser.add_argument("--out_map", help="output module map", type=str, required=True)
    parser.add_argument("--num_modules", help="number of modules to keep", type=int, required=True)
    args = parser.parse_args()
    in_map_fn = args.in_map
    out_map_fn = args.out_map
    max_num_modules = args.num_modules

    in_map_file = codecs.open(in_map_fn, 'r', encoding='utf-8')

    # Read in header
    #
    # modules: 20048
    # modulelinks: 185715
    # nodes: 252373
    # links: 799485
    # codelength: 9.27887
    # *Undirected
    # *Modules 20048
    #
    hdr = {}
    for i in xrange(0, 7):
        ln = in_map_file.readline()
        ln = ln.rstrip()
        f = ln.split()
        if (len(f)==3):
            hdr[f[1]] = f[2]
        elif (len(f)==2):
            hdr[f[0]] = f[1]
        else:
            hdr[f[0]] = ""
    
    # Read in modules
    #
    # 1 "ht-love" 0.188017 0.0581802
    # 
    num_modules = int(hdr['*Modules'])
    print "Found {} modules in file".format(num_modules)
    modules = []
    for i in xrange(1, num_modules+1):
        ln = in_map_file.readline()
        ln = ln.rstrip()
        if (i <= max_num_modules):
            modules.append(ln)

    print "modules: {}".format(modules)

    # Read in nodes
    # 
    # *Nodes 252373
    # 1:1 "ht-love" 0.00618679
    ln = in_map_file.readline()
    if (ln.split()[0]!="*Nodes"):
        print "Format error, expected *Nodes got {}".format(ln)
        exit(1)
    nodes = []
    for ln in in_map_file:
        ln = ln.rstrip()
        mod_num = int(ln.split()[0].split(":")[0])
        if (mod_num <= max_num_modules):
            nodes.append(ln)
        else:
            break

    for ln in in_map_file:
        ln = ln.rstrip()
        if (ln.split()[0]=="*Links"):
            break

    links = []
    for ln in in_map_file:
        ln = ln.rstrip()
        f = ln.split()
        n1 = int(f[0])
        n2 = int(f[1])
        if (n1<=max_num_modules and n2<=max_num_modules):
            links.append(ln)
        
    # Now write it all out
    
    # modules: 20048
    # modulelinks: 185715
    # nodes: 252373
    # links: 799485
    # codelength: 9.27887
    # *Undirected
    # *Modules 20048
    out_map_file = codecs.open(out_map_fn, 'w', encoding='ascii', errors='ignore')
    out_map_file.write('# modules: {}\n'.format(max_num_modules))
    out_map_file.write('# modulelinks: {}\n'.format(len(links)))
    out_map_file.write('# nodes: {}\n'.format(len(nodes)))
    out_map_file.write('# links: {}\n'.format(len(links)))
    out_map_file.write('# codelength: {}\n'.format(hdr['codelength:']))
    out_map_file.write('*Undirected\n')
    out_map_file.write('*Modules {}\n'.format(max_num_modules))
    for ele in modules:
        out_map_file.write("{}\n".format(ele))

    # Nodes
    out_map_file.write("*Nodes {}\n".format(len(nodes)))
    for ele in nodes:
        out_map_file.write(u"{}\n".format(ele).encode('ascii', 'ignore'))

    # Links
    out_map_file.write("*Links {}".format(len(links)))
    for ele in links:
        out_map_file.write("{}\n".format(ele))
    
    out_map_file.close()

