#!/bin/bash
# $1 = input fasta file
# $2 = output alignment file (.aln or .txt)
# $3 = any other options
# example = /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna /data/wrayva/output/extract_regions/DJ/good_dj_alignment.aln
set -e
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK --retree 1 --maxiterate 0 $3 $1 > $2
