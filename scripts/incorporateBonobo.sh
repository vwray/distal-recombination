

module load samtools
module load R

origRegion=DJ
region=DJ
extractDir=/data/wrayva/output/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
sequencesDir=/data/wrayva/output/sequences

extractedRegion=/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa

mkdir -p $regionDir $regionDir/bed $regionDir/sequences $regionDir/paf $regionDir/filtered_paf ${regionDir}/mafft $regionDir/plots/ $regionDir/bonobo/

extractStart=$(cat ${extractDir}/${origRegion}.bed | awk '{print $2}')
extractEnd=$(cat ${extractDir}/${origRegion}.bed | awk '{print $3}')
expectedLength=$((extractEnd-extractStart))

#extract region from chm13 chr22
samtools faidx /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa chr22:${extractStart}-${extractEnd} > $regionDir/chm13chr22${region}.fa

#check for bonobo region B
#sbatch --time=5:00:00 --mem=32g --cpus-per-task=5 ${scriptDir}/run_minimap.sh ${regionDir}/bonobo/bonobo.fna $regionDir/chm13chr22${region}.fa $regionDir/bonobo/bonobo_minimap_${region}_aln.paf
#31751368

sbatch --time=5:00:00 --mem=32g --cpus-per-task=5 ${scriptDir}/run_minimap.sh /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa $regionDir/chm13chr22${region}.fa $regionDir/paf/bonobo_minimap_${region}_aln.paf
#31750974

#filter paf
awk '($17~"tp:A:P")  {print $0}' $regionDir/bonobo/bonobo_minimap_${region}_aln.paf > $regionDir/bonobo/filtered_bonobo_minimap_${region}_aln.paf

#extract bonobo region B to include in tree, add to ${regionDir}/sequences/all_${region}.fna
cat $regionDir/paf/bonobo_minimap_${region}_aln.paf
#chr22:226744-459516	232773	5	232769	+	chr14_pat_hsa13	116126430	8749771	8976845	139590

bonobo_coords=$(cat $regionDir/bonobo/filtered_bonobo_minimap_${region}_aln.paf | sort -nrk10 | head -n1 | awk '{print $6,$8,$9,$5}')

sequenceName=$(echo $bonobo_coords | awk '{print $1}')
echo $sequenceName

bonoboStart=$(echo $bonobo_coords | awk '{print $2}')
echo $bonoboStart

bonoboEnd=$(echo $bonobo_coords | awk '{print $3}')
echo $bonoboEnd

orientation=$(echo $bonobo_coords | awk '{print $4}')
echo $orientation

samtools faidx /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa ${sequenceName}:${bonoboStart}-${bonoboEnd} --reverse-complement --mark-strand no >> ${regionDir}/bonobo/bonobo.fna

samtools faidx ${regionDir}/bonobo/bonobo.fna ${sequenceName}:${bonoboStart}-${bonoboEnd} >> ${regionDir}/good_${region}.fna
#${regionDir}/bonobo/bonobo_${region}.fna

toReplace=$(grep $sequenceName ${regionDir}/bonobo/bonobo.fna)
echo $toReplace

replaceWith=">mPanPan1_$sequenceName"
echo $replaceWith

sed -i "s/$toReplace/$replaceWith/g" ${regionDir}/bonobo/bonobo.fna

#copy bonobo sequence over to the main fasta file
cat ${regionDir}/bonobo/bonobo.fna >> $regionDir/sequences/good_DJ.fna

#run mafft on ${regionDir}/sequences/all_${region}.fna
cd ${regionDir}/mafft
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK --retree 2 --maxiterate 0 ${regionDir}/good_${region}.fna > ${regionDir}/mafft/mafft_${region}_with_bonobo.aln


toReplace=$(grep $sequenceName ${regionDir}/mafft/mafft_${region}.aln)
echo $toReplace

replaceWith=">mPanPan1_$sequenceName"
echo $replaceWith

sed -i "s/$toReplace/$replaceWith/g" ${regionDir}/mafft/mafft_${region}.aln


sbatch --time=5:00:00 --mem=32g --cpus-per-task=20 /data/wrayva/scripts/run_iqtree.sh ${regionDir}/mafft/mafft_${region}.aln ${regionDir}/mafft/ "-o mPanPan1_$sequenceName"
#31806983

#sbatch --time=5:00:00 --mem=32g --cpus-per-task=20 ../regionB2/run_mafft_and_iqtree.sh $region $regionDir mPanPan1_$sequenceName
#31791979

#rename things to have superpop and pop in name
cp ${regionDir}/mafft/mafft_${region}.aln.treefile ${regionDir}/mafft/mafft_${region}_copy.aln.treefile
cd ${genomeDir}
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=${regionDir}/mafft/mafft_${region}_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done

#need to resize bonobo branch
Rscript $scriptDir/plotTrees.r > $regionDir/plots/${region}_big_tree.pdf


#sbatch --time=5:00:00 --mem=32g --cpus-per-task=20 $scriptDir/alignAndCreateTrees50Random.sh
#31727316
#sbatch --time=2:00:00 --mem=16g --cpus-per-task=20 /data/wrayva/scripts/run_iqtree.sh ${subsetDir}/mafft/mafft_${region}.aln ${subsetDir}/mafft/ "-o $(echo $replaceWith | cut -c2-)"
#31731239

#run Parsnp
#sbatch --time=10:00:00 --mem=32g /data/wrayva/scripts/run_parsnp.sh /data/wrayva/output/extract_regions/regionA2/chr22_regionA.fna /data/wrayva/output/extract_regions/regionA2/good_sequences /data/wrayva/output/extract_regions/regionA2/parsnp-output -c
