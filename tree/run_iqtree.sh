#!/bin/bash
# $1 = input alignment file
# $2 = output directory
# $3 = any other options
# example = /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft.fna /data/wrayva/output/extract_regions/regionA/iqtree_output
set -e
module load iqtree
cd $2
iqtree2 -T $SLURM_CPUS_PER_TASK $3 -s $1
