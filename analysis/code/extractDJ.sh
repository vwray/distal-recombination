module load minimap2

cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"
    #get telomere coords from seqtk telo
    teloPosition=$(cat /data/wrayva/output/seqtk_telo/output_rdna_trim_${file%.fna}.txt | head -n1 | awk -F$'\t' '{print $3}')
    #handle empty minimap output
    teloFile="/data/wrayva/output/seqtk_telo/output_rdna_trim_${file%.fna}.txt"; [ -e $teloFile ] && [ -s $teloFile ] || teloPosition=0
    #echo $teloPosition
    echo -e $chrName"\t"0"\t"$teloPosition"\t"telo_${hapName} >> /data/wrayva/output/bed/${genomeName}.bed
    #write telomere coords to bed file

    #run minimap on all sequences with my seq as ref and the DJ as query.
    #minimap2 -x asm20 $file /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa > /data/wrayva/output/minimap_dj/${file%.fna}.paf
    djStart=$(awk '($13~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_dj/${file%.fna}.paf | awk '$12 >= 60'| sort -nrk8 | head -n1 | awk '{print $8}')
    #grab entry 9, tab separated, in the paf file
    djEnd=$(awk '($13~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_dj/${file%.fna}.paf | awk '$12 >= 60'| sort -nrk8 | head -n1 | awk '{print $9}')
    #write DJ coords to bed file
    #echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName}
    echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName} >> /data/wrayva/output/bed/${genomeName}.bed
    # skip trim my sequence at dj end, just record in bed file
    # skip trim my sequence at telomere end, just record in bed file
    # get lengths
    length="$((djEnd - teloPosition))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/bed/lengths.txt
done

#replace sequence identifiers in sequence folder
cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"

    toReplace=$(cat $file | grep -o "\w*sequence_\w*")
    echo $toReplace

    replaceWith=distal_${genomeName}_${chrName}_${hapName}
    echo $replaceWith

    sed -i "s/$toReplace/$replaceWith/g" $file

done


#extract DJ
#module load seqtk
cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"

    djBedLine=$(cat /data/wrayva/output/bed/${genomeName}.bed | grep DJ_${hapName})
    djStart=$(echo $djBedLine | awk '{print $2}')
    djEnd=$(echo $djBedLine | awk '{print $3}')

    #extract DJ
    #write DJ coords to new bedfile
    echo $djBedLine >> /data/wrayva/output/extract_regions/DJ/${file%.fna}.bed

    #extract with seqtk
    #seqtk subseq $file ${file%.fna}.bed > /data/wrayva/output/extract_regions/DJ/dj_${file}

    #extract with cut
    cutStart="$(($djStart - 1))"
    zeroIndexedStart="$(($cutStart - 1))"
    cutEnd="$(($djEnd + 6))"
    zeroIndexedEnd=$cutEnd
    echo ">${file%.fna}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/DJ/dj_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/DJ/dj_${file}

    #output DJ length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/DJ/lengths.txt

done


#put all sequences in one fasta file:
cd /data/wrayva/output/extract_regions/DJ
for file in `ls dj_distal_*fna`; do
    cat $file >> /data/wrayva/output/extract_regions/DJ/all_dj.fna
done

#run muscle
#!/bin/bash
set -e
module load muscle
cd /data/wrayva/output/extract_regions/DJ
muscle -threads $SLURM_CPUS_PER_TASK -align /data/wrayva/output/extract_regions/DJ/all_dj.fna -output /data/wrayva/output/extract_regions/DJ/dj_alignment.fna




sbatch --time=30:00:00 --mem=256g --cpus-per-task=20 /data/wrayva/output/extract_regions/DJ/run_muscle.sh


#try running Muscle on just 50 sequences
cat /data/wrayva/output/extract_regions/DJ/all_dj.fna | head -n100 >> /data/wrayva/output/extract_regions/DJ/first_50_dj.fa

muscle -threads 2 -align /data/wrayva/output/extract_regions/DJ/first_50_dj_chopped.fa -output /data/wrayva/output/extract_regions/DJ/first_50_dj_alignment.fa

#try running mafft-linsi
#!/bin/bash
set -e
module load mafft
mafft-linsi --thread $SLURM_CPUS_PER_TASK /data/wrayva/output/extract_regions/DJ/all_dj.fna > /data/wrayva/output/extract_regions/DJ/mafft-linsi-output.txt

sbatch --time=30:00:00 --mem=128g --cpus-per-task=20 /data/wrayva/output/extract_regions/DJ/run_mafft.sh

#try running parsnp
#!/bin/bash
module load parsnp
parsnp -r /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa -d /data/wrayva/output/extract_regions/DJ/dj_distal_*.fna -o /data/wrayva/output/extract_regions/DJ/parsnp-output/dj_parsnp_output -v

sbatch --time=10:00:00 --mem=32g /data/wrayva/output/extract_regions/DJ/parsnp-output/run_parsnp.sh

sbatch --time=8:00:00 --mem=32g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft.aln /data/wrayva/output/extract_regions/DJ/iqtree_output/


#sanity check dj, checking too short ones
cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 4793795 && $9 >= 4448072' | cut -c-155 | grep distal_HG01993_chr13_haplotype2-0000089

cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | cut -c-155 | grep distal_HG00280_chr15_haplotype1-0000019

djStart=$(awk '($13~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_dj/distal_HG01993_chr13_haplotype2-0000089.paf | awk '$12 >= 60'| sort -nrk8 | sort -nrk10 | head -n1 | awk '{print $8}')
#grab entry 9, tab separated, in the paf file
djEnd=$(awk '($13~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_dj/${file%.fna}.paf | awk '$12 >= 60'| sort -nrk8 | head -n1 | awk '{print $9}')

cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | grep distal_HG00280_chr15_haplotype1-0000019 > /data/wrayva/output/extract_regions/DJ/validation/extracted_dj_distal_HG00280_chr15_haplotype1-0000019.paf

cat /data/wrayva/output/sequences/distal_HG01993_chr13_haplotype2-0000089.fna | cut -c3152804-3152903

#DJ from chr13
#from end
cat /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa | head -n2 | tail -n1 | rev | cut -c-100
cggaagttcaggtcgcagcgttcgtaccgttatggggaagagatgattacaagtttttattcgtcccgtaataccgcgcactgacattagggccggtgag
#from beginning
cat /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa | head -n2 | tail -n1 | cut -c-100
aatgtaaggctgaaagctgtaaaactcctacacaaaaatataagggaaaaatatgcaacgttatgccattgaatttggcagtgggatcttgACTGCCAGC
#DJ from chr22_distal
#from end
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | head -n2 | tail -n1 | rev | cut -c-
CGGAAGTTCAGGTCGCAGCGTTCGTACCGTTATGGGGAAGAGATGATTACAAGTTTTTAATCGTCCCGTAATACCGCGCACTGACATTAGGGCCGGTGAG
#from beginning
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | head -n2 | tail -n1 | cut -c4448074-4448173
AATGTAAGGCTGAAAGCTGTAAAACTCCTACACAAAAATATAAGGGAAAAATATGCAACGTTATGCCATTGAATTTGGCAGTGGGATCTTGACTGCCAGC
#get length of string
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | head -n2 | tail -n1 | tr -d '\n\r' | wc -c
#4793795
#DJ from sample
#from end
cat /data/wrayva/output/extract_regions/DJ/sequences/dj_distal_HG00280_chr15_haplotype1-0000019.fna | head -n2 | tail -n1 | rev | cut -c-100
#from beginning
cat dj_distal_HG00280_chr15_haplotype1-0000019.fna | head -n2 | tail -n1 | cut -c-100

cat /data/wrayva/output/extract_regions/DJ/sequences/dj_distal_HG02004_chr22_haplotype1-0000009.fna | head -n2 | tail -n1 | cut -c-100
HG02004	haplotype1-0000009	chr22

distal_HG01993_chr13_haplotype2-0000089

#create paf with bad djs (the ones that are too short). view in IGV
cat /data/wrayva/output/extract_regions/DJ/too_small_dj.txt | while read -r line; do
    echo $line
    genomeName=$(echo $line | awk '{print $1}')
    chrName=$(echo $line | awk '{print $3}')
    hapName=$(echo $line | awk '{print $2}')
    echo "hapName: $hapName"
    cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | grep distal_${genomeName}_${chrName}_${hapName} >> /data/wrayva/output/extract_regions/DJ/validation/bad_djs.paf
done

#good and bad DJ
for file in `ls /data/wrayva/output/extract_regions/DJ/sequences/dj_distal_*fna`; do
    if [[ $(cat $file | head -n2 | tail -n1 | cut -c-10) == 'AATGTAAGGC' && $(cat $file | head -n2 | tail -n1 | rev | cut -c-10) == 'CGGAAGTTCA' ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/DJ/good_DJ.txt
    else
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/DJ/bad_DJ.txt
    fi
done

#put all good ones in one fasta
for name in `cat /data/wrayva/output/extract_regions/DJ/good_DJ.txt`; do
    cat /data/wrayva/output/extract_regions/DJ/sequences/dj_${name}.fna >> /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna
done

#copy all good ones to separate directory, separate fasta files
for name in `cat /data/wrayva/output/extract_regions/DJ/good_DJ.txt`; do
    cp /data/wrayva/output/extract_regions/DJ/sequences/dj_${name}.fna /data/wrayva/output/extract_regions/DJ/good_sequences/
done

#run muscle on good dj
#!/bin/bash
set -e
module load muscle
cd /data/wrayva/output/extract_regions/DJ
muscle -threads $SLURM_CPUS_PER_TASK -super5 /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna -output /data/wrayva/output/extract_regions/DJ/good_dj_alignment.fna

sbatch --time=30:00:00 --mem=256g --cpus-per-task=20 /data/wrayva/output/extract_regions/DJ/run_muscle_good_dj.sh
#30475115
30480428

#try running parsnp on good dj
#!/bin/bash
module load parsnp
parsnp -r /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa -d /data/wrayva/output/extract_regions/DJ/good_sequences -o /data/wrayva/output/extract_regions/DJ/parsnp-output/good_dj_parsnp_output -v

sbatch --time=10:00:00 --mem=32g /data/wrayva/output/extract_regions/DJ/parsnp-output/run_parsnp_good_dj.sh
#30483682

#make parsnp script general
touch /data/wrayva/scripts/run_parsnp.sh

sbatch --time=10:00:00 --mem=32g /data/wrayva/scripts/run_parsnp.sh /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa /data/wrayva/output/extract_regions/DJ/good_sequences /data/wrayva/output/extract_regions/DJ/parsnp-output

#make muscle script general
touch /data/wrayva/scripts/run_muscle.sh
#!/bin/bash
# $1 = align or super5 (super5 meant for higher mem usage)
# $2 = input fasta file
# $3 = output alignment fasta file
# example = /data/wrayva/scripts/run_muscle.sh align /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna /data/wrayva/output/extract_regions/DJ/good_dj_alignment.fna
set -e
module load muscle
cd /data/wrayva/output/extract_regions/DJ
muscle -threads $SLURM_CPUS_PER_TASK -$1 $2 -output $3

#make mafft script general
touch /data/wrayva/scripts/run_mafft.sh
#!/bin/bash
# $1 = input fasta file
# $2 = output alignment file (.aln or .txt)
# $3 = any other options
# example = /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna /data/wrayva/output/extract_regions/DJ/good_dj_alignment.aln
set -e
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK --retree 1 --maxiterate 0 $3 $1 > $2

sbatch --time=8:00:00 --mem=32g --cpus-per-task=20 /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/DJ/sequences/good_DJ.fna /data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft.aln
#30505126
cat /data/wrayva/output/slurm-30505126.out

#run iqtree
sbatch --time=8:00:00 --mem=32g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft.fna /data/wrayva/output/extract_regions/DJ/iqtree_output
#30508028


#test
module load seqtk
cd /data/wrayva/output/sequences
for file in `ls distal_*fna | head -n30`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"
    djBedLine=$(cat /data/wrayva/output/bed/${genomeName}.bed | grep DJ_${hapName})
    djStart=$(echo $djBedLine | awk '{print $2}')
    djEnd=$(echo $djBedLine | awk '{print $3}')
    echo $djStart
    echo $djEnd



    echo $djBedLine >> /data/wrayva/output/extract_regions/DJ/${file%.fna}.bed

    #extract with seqtk
    #seqtk subseq $file ${file%.fna}.bed > /data/wrayva/output/extract_regions/DJ/dj_${file}
    #extract with cut
    cutStart="$(($djStart - 1))"
    zeroIndexedStart="$(($cutStart - 1))"
    cutEnd="$(($djEnd + 6))"
    zeroIndexedEnd=$cutEnd
    echo ">${chrName}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/DJ/dj_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/DJ/dj_${file}

    #output DJ length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/DJ/lengths.txt
done


cat /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa | head -n2 | tail -n1 | rev | cut -c-100
cat dj_distal_HG00099_chr14_haplotype2-0000155.fna | head -n2 | tail -n1 | rev | cut -c-100

#tail end missing 6 characters

cat /data/wrayva/output/bed/lengths.txt | grep chr21 | awk '$4 >= 2714389 && $4 <= 3312526' | wc -l


file=distal_NA20805_chr13_haplotype2-0000136
distal_NA20805_chr13_haplotype2-0000136
NA20805	haplotype2-0000136	chr13	32673
 - 3281

file=distal_HG02280_chr15_haplotype2-0000203
chr15	3436301	3497864	DJ_haplotype2-0000203
file=NA20905_chr22_haplotype2-0000055.fna



#chr22 	 2701070 3046849
#chr22	2701070 3046849
#chr22	2701070	3046849





Tree file: good_dj_alignment_mafft_copy.aln.treefile

#super pop labels
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')
    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}
    sed -i "s/$toReplace/$replaceWith/g" $file
done

#color by super population
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy.aln.treefile")
pops <- list(EUR=c(tree$tip.label[grep("EUR", tree$tip.label)]),
             EAS=c(tree$tip.label[grep("EAS", tree$tip.label)]),
             AMR=c(tree$tip.label[grep("AMR", tree$tip.label)]),
             AFR=c(tree$tip.label[grep("AFR", tree$tip.label)]),
             SAS=c(tree$tip.label[grep("SAS", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, pops,
                       group_name = "population")
ggtree(grouped_tree2, aes(color=population), layout='circular') +
  theme(legend.position=c(.95,.95), legend.key.size = unit(.2, 'cm')) +
  geom_tiplab(align=TRUE, linesize=.1, size=1)


chrs <- list(chr13=c(tree$tip.label[grep("chr13", tree$tip.label)]),
           chr14=c(tree$tip.label[grep("chr14", tree$tip.label)]),
           chr15=c(tree$tip.label[grep("chr15", tree$tip.label)]),
           chr21=c(tree$tip.label[grep("chr21", tree$tip.label)]),
           chr22=c(tree$tip.label[grep("chr22", tree$tip.label)])
)
grouped_tree3 <- groupOTU(tree, chrs,
                     group_name = "chromosome")
ggtree(grouped_tree3, aes(color=chromosome), layout='circular') +
  theme(legend.position=c(.95,.95), legend.key.size = unit(.1, 'cm')) +
  geom_tiplab(align=TRUE, linesize=.1, size=1)


ggsave("/data/wrayva/output/plots/djTree_circular_color_by_super_pop_and_chr.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)



#get number of samples from each super population
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    #file=/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy.aln.treefile
    #toReplace=distal_$genomeName

    ls /data/wrayva/output/sequences/ | grep

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')
    echo $superpop

    #replaceWith=${superpop}_${pop}_${genomeName}
    #sed -i "s/$toReplace/$replaceWith/g" $file
done

cd /data/wrayva/output/sequences
cd /data/wrayva/output/extract_regions/regionA2
eur_count=0
eas_count=0
amr_count=0
afr_count=0
sas_count=0
#for file in `ls distal_*fna`; do
for file in `cat good_regionA.txt`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')
    echo $superpop
    if [[ $superpop == "EUR" ]]; then
        ((eur_count++))
    fi
    if [[ $superpop == "EAS" ]]; then
        ((eas_count++))
    fi
    if [[ $superpop == "AMR" ]]; then
        ((amr_count++))
    fi
    if [[ $superpop == "AFR" ]]; then
        ((afr_count++))
    fi
    if [[ $superpop == "SAS" ]]; then
        ((sas_count++))
    fi
done
echo EUR:$eur_count EAS:$eas_count AMR:$amr_count AFR:$afr_count SAS:$sas_count



Total: 391 haplotypes  EUR:31 EAS:61 AMR:113 AFR:106 SAS:80
DJ: 267 haplotypes  EUR:25 EAS:35 AMR:63 AFR:79 SAS:65
Region A: 272 haplotypes  EUR:19 EAS:42 AMR:83 AFR:68 SAS:60


The Distal Regions of Human Acrocentric Short Arms are not Chromosome Specific





#get average number of mutations per column of multi-alignment
#convert MSA to format with name space aligned seq on each line
cd /data/wrayva/output/extract_regions/DJ
cat good_dj_alignment_mafft_copy.aln

file=good_dj_alignment_mafft_copy.aln
toReplace=$(cat $file | grep -o ":\w*-\w*")
echo $toReplace

replaceWith=' '
echo $replaceWith

for toReplace2 in `cat $file | grep -o ":\w*-\w*"`; do
    sed -i "s/$toReplace2/$replaceWith/g" $file
done

cat $file | grep ">*:"

for seqIdent in `cat $file | grep ">*:"`; do
    beg=$(echo $seqIdent | cut -d ":" -f1)
    sed -i "s/$seqIdent/$beg/g" $file
done



#align DJ to bonobo haplotype

#!/bin/bash
set -e
module load minimap2
minimap2 -cx asm20 /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa > /data/wrayva/output/extract_regions/DJ/bonobo/minimap_chr4_regA_bonobo.paf

sbatch --time=5:00:00 --mem=24g --cpus-per-task=5 /data/wrayva/output/extract_regions/DJ/bonobo/run_minimap_dj_bonobo.sh
#31416638

#filter paf
cd /data/wrayva/output/extract_regions/DJ/bonobo
awk '($17~"tp:A:P")  {print $0}' minimap_chr4_dj_bonobo.paf > filtered_minimap_chr4_dj_bonobo.paf

cat filtered_minimap_chr4_dj_bonobo.paf | sort -nrk10 | head -n4 | cut -c-150

#extract bonobo dj from minimap alignment coords
coords=$(cat filtered_minimap_chr4_dj_bonobo.paf | sort -nrk10 | head -n1 | awk '{print $1,$3,$4,$5}')

#extracting the first ~20kb of the DJ sequences to see if there are a lot of identical sequences on that region
module load samtools
cd /data/wrayva/output/extract_regions/DJ
for file in `grep ">" good_dj_alignment_mafft.aln`; do
    echo ${file#>}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"
    #echo $genomeName
    #echo $chrName
    echo $file
    echo ${file#>}:1-20000

    samtools faidx good_dj_alignment_mafft.aln ${file#>}:1-20000 >> good_dj_alignment_mafft_copy_first20k2.aln
done

sbatch --time=5:00:00 --mem=8g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy_first20k.aln /data/wrayva/output/extract_regions/DJ/
#31427184

good_dj_alignment_mafft_copy_first20k_copy.aln.treefile

#rename to include superpop and pop
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy_first20k_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done
cd /data/wrayva/output/extract_regions/DJ

#create tree plots
module load R
R

library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy_first20k_copy.aln.treefile")
chrs <- list(chr13=c(tree$tip.label[grep("chr13", tree$tip.label)]),
             chr14=c(tree$tip.label[grep("chr14", tree$tip.label)]),
             chr15=c(tree$tip.label[grep("chr15", tree$tip.label)]),
             chr21=c(tree$tip.label[grep("chr21", tree$tip.label)]),
             chr22=c(tree$tip.label[grep("chr22", tree$tip.label)]),
             chr4=c(tree$tip.label[grep("chr4", tree$tip.label)]),
             bonobo=c(tree$tip.label[grep("chr3", tree$tip.label)])
)
grouped_tree1 <- groupOTU(tree, chrs,
                       group_name = "chromosome")
ggtree(grouped_tree1, aes(color=chromosome), size=1, layout='circular') +
  theme(legend.position=c(.97,.97), legend.text = element_text(size = 10), legend.key.size = unit(.5, 'cm'),plot.margin=margin(200, 200, 200, 200)) +
  geom_tiplab(align=TRUE, linesize=.5, size=2, fontface='bold')

ggsave("/data/wrayva/output/plots/djTree_wfirst20k_color_by_chr.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

pops <- list(EUR=c(tree$tip.label[grep("EUR", tree$tip.label)]),
             EAS=c(tree$tip.label[grep("EAS", tree$tip.label)]),
             AMR=c(tree$tip.label[grep("AMR", tree$tip.label)]),
             AFR=c(tree$tip.label[grep("AFR", tree$tip.label)]),
             SAS=c(tree$tip.label[grep("SAS", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, pops,
                       group_name = "population")
ggtree(grouped_tree2, aes(color=population), size=1, layout='circular') +
  theme(legend.position=c(.97,.97), legend.text = element_text(size = 10), legend.key.size = unit(.5, 'cm'),plot.margin=margin(200, 200, 200, 200)) +
  geom_tiplab(align=TRUE, linesize=.5, size=2, fontface='bold')
ggsave("/data/wrayva/output/plots/djTree_wfirst20k_color_by_super_pop.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

#get first 20kb of all sequences
for name in `cat /data/wrayva/output/extract_regions/DJ/good_DJ.txt`; do
    cp /data/wrayva/output/extract_regions/DJ/sequences/dj_${name}.fna /data/wrayva/output/extract_regions/DJ/good_sequences/
done

cd /data/wrayva/output/extract_regions/DJ/good_sequences/
for file in `ls dj_distal_*`; do
    file2=${file#dj_distal_}
    #file2=${file2%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f3)
    chrName=$(echo ${file%.fna} | cut -d '_' -f4)
    hapName=$(echo ${file%.fna} | cut -d '_' -f5)
    echo "hapName: $hapName"
    echo $genomeName
    echo $chrName
    echo $file
    echo $file2

    line=$(cat $file | head -n1 | cut -d '>' -f2)

    samtools faidx $file ${line}:1-20000 > /data/wrayva/output/extract_regions/DJ/sequences_first20kb/$file2
done

#extract first 20kb of ref
samtools faidx /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa chr13:5424523-5770548:1-20000 > /data/wrayva/output/extract_regions/DJ/ref_first20kb.fa

#run parsnp on first 20kb with -c option to not throw away any sequences
sbatch --time=10:00:00 --mem=32g /data/wrayva/scripts/run_parsnp.sh /data/wrayva/output/extract_regions/DJ/ref_first20kb.fa /data/wrayva/output/extract_regions/DJ/sequences_first20kb /data/wrayva/output/extract_regions/DJ/parsnp-first20kb-output -c
#31458168


#convert MSA to VCF to run ARG method
source myconda
conda activate dfv
msa2vcf good_dj_alignment_mafft_copy.aln distal_HG00099_chr13_haplotype1-0000019:2235491-2579860

#above did not work the way I want it to. trying the below
snp-sites -v -o good_dj_alignment_mafft_copy_vcf.vcf good_dj_alignment_mafft_copy.aln

#convert MSA with just 50 samples
module load snp-sites
cd /data/wrayva/output/extract_regions/DJ/50random/mafft
snp-sites -v -o mafft_DJ_vcf.vcf mafft_DJ.aln

cat /data/wrayva/output/extract_regions/DJ/50random/mafft/mafft_DJ_vcf.vcf | head -n10
cat /data/wrayva/output/extract_regions/DJ/50random/mafft/mafft_DJ.aln | head -n2 | tail -n1 | cut -c-500


#!/bin/bash
set -e
numberOfSamples=51
diploidSamples=51
sequenceLength="366_818"
sequenceLenInt=366818
recombRate=1e-8
mutationRate=1.2e-8
baseDirectory="/data/wrayva/gitRepos/args"
scriptRepo="${baseDirectory}/scripts"
rentPlusDir=/data/wrayva/gitRepos/RentPlus

outputDirectory=/data/wrayva/output/extract_regions/DJ/50random/args/2024-08-06T13:38:39_366_818

vcf_file=/data/wrayva/output/extract_regions/DJ/50random/mafft/mafft_DJ_vcf.vcf

#run RENT+ with branch lengths
cd ${outputDirectory}
module load java
java -Xmx24g -jar $rentPlusDir/RentPlus.jar -t rentOutput.txt

cd /data/wrayva/output/extract_regions/DJ/50random/args
sbatch --time=8:00:00 --mem=32g --cpus-per-task=20 run_rentPlus.sh
#32331350


Rscript $scriptDir/plotSmallTrees.r rentTree10kb.txt rentTree10kb.pdf rentTree10kb2.pdf
