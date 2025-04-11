#!/bin/bash
set -e

#run with
#sbatch --time=8:00:00 --mem=32g --cpus-per-task=20 $scriptDir/run_iqtree_sections.sh

origRegion=DJ
region=DJ
prefix=dj
extractDir=/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts
#/data/wrayva/scripts
#genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
#sequencesDir=/data/wrayva/output/sequences
#subsetDir=$regionDir/50random

cd ${extractDir}/DJ

module load iqtree
for file in `cat rfdist3/msas.txt`; do
    range2=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
    echo $file
    echo $range2
    iqtree2 -T $SLURM_CPUS_PER_TASK -s rfdist3/${file}
# removing outgroup for now
# -o mPanPan1_chr14_pat_hsa13_${range2}
done