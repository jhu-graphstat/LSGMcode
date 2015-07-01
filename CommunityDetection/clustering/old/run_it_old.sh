#!/bin/sh

scripts/gdf_to_edge.py

# Inputs
edge_list=graphs/sn_eg_edge.txt
sigma=0.5
lambda=1.0 # smaller lambda -- less clusters, larger lambda -- more clusters
thresh=0 # admit all
trials=10

# Setup
gawk -F',' '{print $1; print $2}' $edge_list | sort -u > tmp/v.txt
scripts/edge_to_infomap.pl $edge_list tmp/v.txt $sigma $thresh > tmp/out.net

# Run infomap
# scripts/infomap 345234 $NET $TRIALS | grep "Done"
scripts/infomap 345234 tmp/out.net $trials $lambda 

tail -n +2 tmp/out.clu | paste -d" " tmp/v.txt - > tmp/clusters.txt

scripts/label_clusters.py 
