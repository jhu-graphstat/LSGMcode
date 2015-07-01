#!/bin/sh

in_graph=tmp_workshop/twitter_sg0.gpckl

echo Writing weighted graph ...
scripts/write_weighted_graph.py --in_graph $in_graph --out_graph tmp/g.net --out_label tmp/g_lbl.txt

echo Performing infomap ...
scripts/infomap 1234 tmp/g.net 2 1.0
