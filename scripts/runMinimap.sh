#!/bin/bash

module load minimap2
file='/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample/HG00097/verkko-hi-c/analysis/HG00097.refOriented.fasta.gz'
minimap2 -cx asm20 $file /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa > /data/wrayva/output/rerun_v4/HG00097_DJ.paf
