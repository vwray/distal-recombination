#!/bin/bash
set -e

region=regionA
extractDir=/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts
#/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
#sequencesDir=/data/wrayva/output/sequences
#subsetDir=$regionDir/50random
msaDir=${regionDir}/chunk
msaFile="/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/regionA/good_regionA_alignment_mafft_retree2.aln"

#chunk MSA
module load samtools
cd ${regionDir}
chunkSize=5000
windowSize=2500

#get MSA length
#awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $msaFile


msalength=`awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;} { seqlen += length($0)}END{print seqlen}' $msaFile | head -n1`



#msalength=366783
i=1
while [ $((${i}+${chunkSize})) -le $msalength ]
do
    for file in `grep ">" $msaFile`; do
    #for file in `grep ">" good_dj_alignment_mafft.aln | grep -v "distal_HG01981_chr21_haplotype1" | grep -v "distal_HG00621_chr14_haplotype1"`; do
        #echo ${file#>}
        #genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
        #chrName=$(echo ${file%.fna} | cut -d '_' -f3)
        #hapName=$(echo ${file%.fna} | cut -d '_' -f4)
        #echo "hapName: $hapName"
        #echo $genomeName
        #echo $chrName
        #echo $file
        range=${i}-$((${i}+${chunkSize}-1))
        echo ${file#>}:$range

        samtools faidx $msaFile ${file#>}:$range >> ${msaDir}/aln/aln_${range}.aln

        filename=aln_${range}.aln

    done
    echo $filename >> ${msaDir}/msas.txt
    i=$((i+${windowSize}))
done

#testing for one file, one sequence
#file=`grep ">" good_dj_alignment_mafft.aln | head -n1`
#range=${i}-$((${i}+${windowSize}-1))
#echo ${file#>}:$range
#samtools faidx good_dj_alignment_mafft_copy.aln ${file#>}:$range >> rfdist_test/dj_aln_${range}.aln
#samtools faidx good_dj_alignment_mafft.aln ${file#>}:$range >> rfdist_test/non_copy_dj_aln_${range}.aln

#test
#grep ">" good_dj_alignment_mafft.aln | grep -v "distal_HG01981_chr21_haplotype1" | grep -v "distal_HG00621_chr14_haplotype1"
#267
#exclude sequences 	distal_HG01981_chr21_haplotype1 and distal_HG00621_chr14_haplotype1


#run iqtree on regions
#module load iqtree
#for file in `ls rfdist`; do
#    range2=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
#    echo $file
#    echo $range2
#    iqtree2 -T $SLURM_CPUS_PER_TASK -o mPanPan1_chr14_pat_hsa13_${range2} -s rfdist/${file}
#done

#add aln files to a txt file
#cd $msaDir
#chunkSize=5000
#windowSize=2500
#msalength=366783
#i=1
#while [ $((${i}+${chunkSize})) -le $msalength ]
#do
#    range=${i}-$((${i}+${chunkSize}-1))
#    filename=dj_aln_${range}.aln
#    echo $filename >> msas.txt
#    i=$((i+${windowSize}))
#done

sbatch --time=8:00:00 --mem=32g --cpus-per-task=20 $scriptDir/run_iqtree_sections.sh
#39872157

#submit slurm of jobs
cd $msaDir
touch file.swarm
for file in `cat msas.txt`; do
    range2=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
    echo $file
    echo $range2
    echo "cd $msaDir/aln ; iqtree2 -T \$SLURM_CPUS_PER_TASK -s ${msaDir}/aln/${file}" >> $msaDir/file.swarm
    #iqtree2 -T $SLURM_CPUS_PER_TASK -s rfdist3/${file}
done

swarm --module iqtree -g1 file.swarm
#39944859
40702369

#process trees
#rename things to have superpop and pop in name
#chunkSize=5000
#windowSize=2500
#msalength=366783

mkdir ${msaDir}/treefilecopies
echo "Renaming haplotypes to have superpop in name"
cd ${genomeDir}
i=1
while [ $((${i}+${chunkSize})) -le $msalength ]
do
    range=${i}-$((${i}+${chunkSize}-1))
    echo $range

    cp ${msaDir}/aln/aln_${range}.aln.treefile ${msaDir}/treefilecopies/aln_${range}_copy.aln.treefile
    file=${msaDir}/treefilecopies/aln_${range}_copy.aln.treefile

    for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
        toReplace=distal_${genomeName}

        pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
        superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

        replaceWith=${superpop}_${pop}_${genomeName}

        sed -i "s/$toReplace/$replaceWith/g" $file
    done

    #also replace haplotype2-0000156_2672465-3016694_105001-125000: with hap2:
    hap1ReplaceWith=hap1:
    hap2ReplaceWith=hap2:

    #need to do this inside a for loop; grep returns multiple
    for hap1ToReplace in `cat $file | grep -o "haplotype1-[0-9]*_[0-9]*-[0-9]*_[0-9]*-[0-9]*:"`; do
        sed -i "s/$hap1ToReplace/$hap1ReplaceWith/g" $file
    done
    for hap2ToReplace in `cat $file | grep -o "haplotype2-[0-9]*_[0-9]*-[0-9]*_[0-9]*-[0-9]*:"`; do
        sed -i "s/$hap2ToReplace/$hap2ReplaceWith/g" $file
    done


    #hap1ToReplace=`cat $file | grep "haplotype1-[0-9]*_[0-9]*-[0-9]*_[0-9]*-[0-9]*:"`
    #hap2ToReplace=`cat $file | grep "haplotype2-[0-9]*_[0-9]*-[0-9]*_[0-9]*-[0-9]*:"`
    #sed -i "s/$hap1ToReplace/$hap1ReplaceWith/g" $file
    #sed -i "s/$hap2ToReplace/$hap2ReplaceWith/g" $file

    i=$((i+${windowSize}))
done

#testing out with one tree file
#file=/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/DJ/rfdist2/dj_aln_1-20000.aln.treefile
#hapToReplace=`cat $file | grep "haplotype2-[0-9]+_[0-9]+-[0-9]+_[0-9]+-[0-9]+:"`
#hap2ToReplace=`cat $file | grep "haplotype2-[0-9]*_[0-9]*-[0-9]*_[0-9]*-[0-9]*:"`

#put tree file names in file
cd $msaDir/treefilecopies
echo "TreeFileNames" >> $msaDir/treefilenames.csv
echo "Positions" >> $msaDir/treefilepositions.csv
for file in `ls | sort -n -k1.5`; do
    echo $file >> $msaDir/treefilenames.csv
    position=`echo $file | cut -d '_' -f2 | cut -d '-' -f2`
    newPosition=$(( position - 2500 ))
    echo $newPosition >> $msaDir/treefilepositions.csv
done

cd $msaDir
n=1
while [ $n -le 141 ]
do
    if [ $((n % 2)) -eq 1 ]; then
        echo $n
        cat treefilepositions.csv | head -$n | tail -1 >> treefilepositions_odd.csv
        cat treefilenames.csv | head -$n | tail -1 >> treefilenames_odd.csv
    fi
    n=$((n+1))
done


#testing
#ls | sort -n -k1.8

#cd $msaDir
#skip the first line, which is the header
#numlines=`wc -l treefilenames.csv | cut -d ' ' -f1`
#for file in `cat treefilenames.csv | tail -n$(( numlines - 1 ))`; do
    #get the position
#    echo $file | cut -d '_' -f3 | cut -d '-' -f2
#done


#compute RF distance between trees
cd $msaDir
computeRFDistance.r
