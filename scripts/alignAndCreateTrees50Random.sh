#!/bin/bash
set -e

#run with
#sbatch --time=8:00:00 --mem=32g --cpus-per-task=20 $scriptDir/alignAndCreateTrees50Random2.sh
#31911854

module load samtools
module load R

origRegion=DJ
region=DJ
prefix=dj
extractDir=/data/wrayva/output/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
sequencesDir=/data/wrayva/output/sequences
subsetDir=$regionDir/50random

cd $regionDir
mkdir -p $subsetDir ${subsetDir}/mafft

#for each good haplotype, get its superpopulation and add haplotype to the superpopulation file
echo "Sorting haplotypes by their superpopulation"
cd ${regionDir}
for file in `cat good_${region}.txt`; do
    genomeName=$(echo ${file} | cut -d '_' -f2)
    chrName=$(echo ${file} | cut -d '_' -f3)
    hapName=$(echo ${file} | cut -d '_' -f4)

    pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $file >> ${subsetDir}/${superpop}_haps.txt
done

echo "Choosing 10 random haplotypes from each superpopulation"
cd ${subsetDir}
for hapFile in `ls *_haps.txt`; do
    for file in `sort -R $hapFile | head -n10`; do
        echo $file >> ten_haps_from_each_pop.txt
    done
done
echo "Collecting 10 random haplotypes into one fasta file"
#put all random samples in one fasta
for name in `cat ten_haps_from_each_pop.txt`; do
    cat ${regionDir}/sequences/${prefix}_${name}.fna >> ${subsetDir}/random_samples.fna
done

#echo "Retrieving bonobo sample"
#bonobo_coords=$(cat $regionDir/paf/filtered_bonobo_minimap_${region}_aln.paf | sort -nrk10 | head -n1 | awk '{print $6,$8,$9,$5}')

#sequenceName=$(echo $bonobo_coords | awk '{print $1}')
#bonoboStart=$(echo $bonobo_coords | awk '{print $2}')
#bonoboEnd=$(echo $bonobo_coords | awk '{print $3}')
#orientation=$(echo $bonobo_coords | awk '{print $4}')

#samtools faidx /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa ${sequenceName}:${bonoboStart}-${bonoboEnd} >> ${subsetDir}/random_samples.fna

cat ${regionDir}/bonobo/bonobo.fna >> ${subsetDir}/random_samples.fna

echo "Running mafft"
#run mafft on random samples
cd ${subsetDir}/mafft
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK --retree 2 --maxiterate 0 ${subsetDir}/random_samples.fna > ${subsetDir}/mafft/mafft_${region}.aln
#--auto

#toReplace=$(grep $sequenceName ${subsetDir}/random_samples.fna)
#replaceWith=">mPanPan1_$sequenceName"

#sed -i "s/$toReplace/$replaceWith/g" ${subsetDir}/mafft/mafft_${region}.aln

echo "Running IQ-TREE"
module load iqtree
/data/wrayva/scripts/run_iqtree.sh ${subsetDir}/mafft/mafft_${region}.aln ${subsetDir}/mafft/ "-o $(echo $replaceWith | cut -c2-)"
iqtree2 -T $SLURM_CPUS_PER_TASK -o $(grep "mPanPan" ${subsetDir}/mafft/mafft_${region}.aln | cut -c2-) -s ${subsetDir}/mafft/mafft_${region}.aln

echo "Renaming haplotypes to have superpop in name"
#rename things to have superpop and pop in name
cp ${subsetDir}/mafft/mafft_${region}.aln.treefile ${subsetDir}/mafft/mafft_${region}_copy.aln.treefile
cd ${genomeDir}
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=${subsetDir}/mafft/mafft_${region}_copy.aln.treefile
    toReplace=${prefix}_distal_$genomeName

    pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done

#need to resize bonobo branch
#cat ${regionDir}/good_dj_alignment_mafft.aln.treefile | grep "Pan"

echo "Plotting trees"
Rscript $scriptDir/plotTrees.r ${subsetDir}/mafft/mafft_${region}_1-10000.aln.treefile $regionDir/plots/${region}_small_tree.pdf




#chop MSA into segments
module load samtools
for file in `grep ">" ${subsetDir}/mafft/mafft_${region}.aln`; do
    echo ${file#>}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f1 | cut -d '>' -f2) #orig file had -f3
    hapName=$(echo ${file%.fna} | cut -d '_' -f3) #orig file had -f4
    #echo "hapName: $hapName"
    #echo $genomeName
    #echo $chrName

    samtools faidx ${subsetDir}/mafft/mafft_${region}.aln ${file#>}:1-10000 >> ${subsetDir}/mafft/mafft_${region}_1-10000.aln
    samtools faidx ${subsetDir}/mafft/mafft_${region}.aln ${file#>}:10000-20000 >> ${subsetDir}/mafft/mafft_${region}_10000-20000.aln
    samtools faidx ${subsetDir}/mafft/mafft_${region}.aln ${file#>}:20000-30000 >> ${subsetDir}/mafft/mafft_${region}_20000-30000.aln
    samtools faidx ${subsetDir}/mafft/mafft_${region}.aln ${file#>}:30000-40000 >> ${subsetDir}/mafft/mafft_${region}_30000-40000.aln
    samtools faidx ${subsetDir}/mafft/mafft_${region}.aln ${file#>}:40000-50000 >> ${subsetDir}/mafft/mafft_${region}_40000-50000.aln
done

#outgroup: mPanPan1_chr14_pat_hsa13

#!/bin/bash
set -e

origRegion=DJ
region=DJ
prefix=dj
extractDir=/data/wrayva/output/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
sequencesDir=/data/wrayva/output/sequences
subsetDir=$regionDir/50random

cd $regionDir
module load iqtree
iqtree2 -T $SLURM_CPUS_PER_TASK -o mPanPan1_chr14_pat_hsa13_1-10000 -s ${subsetDir}/mafft/mafft_${region}_1-10000.aln
iqtree2 -T $SLURM_CPUS_PER_TASK -o mPanPan1_chr14_pat_hsa13_10000-20000 -s ${subsetDir}/mafft/mafft_${region}_10000-20000.aln
iqtree2 -T $SLURM_CPUS_PER_TASK -o mPanPan1_chr14_pat_hsa13_20000-30000 -s ${subsetDir}/mafft/mafft_${region}_20000-30000.aln
iqtree2 -T $SLURM_CPUS_PER_TASK -o mPanPan1_chr14_pat_hsa13_30000-40000 -s ${subsetDir}/mafft/mafft_${region}_30000-40000.aln
iqtree2 -T $SLURM_CPUS_PER_TASK -o mPanPan1_chr14_pat_hsa13_40000-50000 -s ${subsetDir}/mafft/mafft_${region}_40000-50000.aln


cd ${subsetDir}/mafft
for file in `ls mafft_*.treefile`; do
    cp $file ${file%.aln.treefile}_copy.aln.treefile
done
cp ${regionDir}/mafft/mafft_${region}.aln.treefile ${regionDir}/mafft/mafft_${region}_copy.aln.treefile


file1=${subsetDir}/mafft/mafft_DJ_1-10000_copy.aln.treefile
file2=${subsetDir}/mafft/mafft_DJ_10000-20000_copy.aln.treefile
file3=${subsetDir}/mafft/mafft_DJ_20000-30000_copy.aln.treefile
file4=${subsetDir}/mafft/mafft_DJ_30000-40000_copy.aln.treefile
file5=${subsetDir}/mafft/mafft_DJ_40000-50000_copy.aln.treefile
cd ${genomeDir}
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    toReplace=distal_$genomeName

    pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file1
    sed -i "s/$toReplace/$replaceWith/g" $file2
    sed -i "s/$toReplace/$replaceWith/g" $file3
    sed -i "s/$toReplace/$replaceWith/g" $file4
    sed -i "s/$toReplace/$replaceWith/g" $file5
done

echo $file1 | grep

#sed -i "s/:/_/g" ${regionDir}/run_iqtree_sections.sh
cd ${subsetDir}/mafft
for file in `ls mafft_DJ_*_copy.aln.treefile`; do
    echo $file
    positions=$(echo $file | grep -o '[0-9]*-[0-9]*')
    toReplace=$(grep -o "mPanPan1_chr14_pat_hsa13_${positions}:0.[0-9]*" $file)
    replaceWith="mPanPan1_chr14_pat_hsa13_${positions}:0.0000000000001"
    sed -i "s/$toReplace/$replaceWith/g" $file

    for hapString in `grep -o 'haplotype[1-2]-[0-9]\+_[0-9]\+-[0-9]\+' $file`; do
        toReplace=$hapString
        echo $toReplace
        replaceWith=hap$(echo $hapString | cut -c 10)
        echo $replaceWith
        sed -i "s/$toReplace/$replaceWith/g" $file
    done
    sed -i "s/_${positions}//g" $file
done

cd ${subsetDir}/mafft
for file in `ls mafft_DJ_*_copy.aln.treefile`; do
    Rscript $scriptDir/plotSmallTrees.r $file $regionDir/plots/${file%.aln.treefile}_by_chr4.pdf $regionDir/plots/${file%.aln.treefile}_by_pop4.pdf
done
