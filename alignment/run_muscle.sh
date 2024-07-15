#!/bin/bash
# $1 = align or super5 (super5 meant for higher mem usage)
# $2 = input fasta file
# $3 = output alignment fasta file
# example = /data/wrayva/scripts/run_muscle.sh align /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna /data/wrayva/output/extract_regions/DJ/good_dj_alignment.fna
set -e
module load muscle
cd /data/wrayva/output/extract_regions/DJ
muscle -threads $SLURM_CPUS_PER_TASK -$1 $2 -output $3
