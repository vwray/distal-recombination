#!/bin/bash
# $1 = ref fasta
# $2 = input folder containing separate fastas
# $3 = output folder
#example = /data/wrayva/scripts/run_parsnp.sh /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa /data/wrayva/output/extract_regions/DJ/good_sequences /data/wrayva/output/extract_regions/DJ/parsnp-output/good_dj_parsnp_output
module load parsnp
parsnp -r $1 -d $2 -o $3 -v
