#!/usr/bin/env bash
realn='./testdat/realn.bam'
minibam='./testdat/miniBAM.bam'
U1='./testdat/U1U11.bed'
out='./testout/test'

## Prepare for realignment
## Run realignment
## Post-process realignments
python ./src/postprocess_realign.py $realn $minibam $out $U1
## Calculate weighted depth for EB filter
python ./src/calc_weighted_depth.py ${out}_processed.bam ${out}_mstat.csv ${out}_umap.txt $U1 test_donor ${out}_eb_depth.tsv