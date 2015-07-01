import networkx as nx
import numpy as np
import scipy.io as sio
import csv

# construct the first graph    
G1 = nx.DiGraph()
efile = open('a_w50000_edges.txt','r')
for line in efile:
    tmp = [word for word in line.split()]
    if len(tmp) == 5 and tmp[2]>0: 
        #G1.add_edge(int(tmp[0]),int(tmp[1]),weight=float(tmp[2]))
        G1.add_edge(int(tmp[0]),int(tmp[1]),weight=1.0)

# construct the second graph
G2 = nx.DiGraph()
efile = open('b_w50000_edges.txt','r')
for line in efile:
    tmp = [word for word in line.split()]
    if len(tmp) == 5 and tmp[2]>0: 
        #G2.add_edge(int(tmp[0]),int(tmp[1]),weight=float(tmp[2]))
        G2.add_edge(int(tmp[0]),int(tmp[1]),weight=1.0)

# read the first file
file = open('a_w50000_nodes.txt','r')
dict1 = {}
for line in file:
    tmp = [word for word in line.split()]
    dict1[tmp[1]]=int(tmp[0])      

# read the second file
file = open('b_w50000_nodes.txt','r')
dict2 = {}
for line in file:
    tmp = [word for word in line.split()]
    dict2[tmp[1]]=int(tmp[0])      

# find matching nodes
a = set(dict1.keys())
b = set(dict2.keys())
sharedNodes = a.intersection(b) # this list is bigger than the list of nodes in the graphs
match_a = set([dict1[n] for n in list(sharedNodes)])
a = match_a - match_a.intersection(set(G1.nodes())) # find the nodes that are not in G1
a =  [dict1.keys()[dict1.values().index(i)] for i in list(a)]
match_b = set([dict2[n] for n in list(sharedNodes)])
b = match_b - match_b.intersection(set(G2.nodes()))
b = [dict2.keys()[dict2.values().index(i)] for i in list(b)] # find ithe nodes that are not in G2
sharedNodes = sharedNodes - set(a) - set(b)
match_a = [dict1[n] for n in list(sharedNodes)]
match_b = [dict2[n] for n in list(sharedNodes)]
H1 = G1.subgraph(match_a)
H2 = G2.subgraph(match_b)
A = nx.adjacency_matrix(H1,match_a)
B = nx.adjacency_matrix(H2,match_b)
sio.savemat('A_binary.mat',{'A':A})
sio.savemat('B_binary.mat',{'B':B})
