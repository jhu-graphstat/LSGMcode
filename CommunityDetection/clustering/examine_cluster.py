#!/usr/bin/env python

import codecs
import sys

lbl_fn = 'results/v1/g_lbl.txt'
clust_fn = 'results/v1/g.tree'
top_num = 10

# Create node number to label dictionary
lbl_file = codecs.open(lbl_fn, 'r', encoding='utf-8')
lbl_dict = {}
for ln in lbl_file:
    ln = ln.rstrip()
    f = ln.split()
    lbl_dict[int(f[0])] = f[1]
lbl_file.close()

# Now get desired cluster
cnum = 0
if (len(sys.argv)>1):
    cnum = int(sys.argv[1])
print "Info for cluster: {}".format(cnum)

# Read in cluster file
count = 0
clust_file = open(clust_fn, 'r')
for ln in clust_file:
    ln = ln.rstrip()
    if (ln[0]=='#'):
        continue
    f = ln.split()
    cc = int(f[0].split(':')[0])
    if (cc==cnum):
        if (count < top_num):
            nnum = int(f[2][6:-1])
            nname = lbl_dict[nnum]
            print 'Line: {}, node: {}'.format(ln, nname)
            count += 1
        else:
            break
clust_file.close()


