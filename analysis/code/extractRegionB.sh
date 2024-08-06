genomeName=HG00099

#test bad regionB
file=distal_HG00099_chr14_haplotype1-0000022.fna

ct 2-5: 232772
ct 3-5: 229363
ct 4-5: 207637
reg C = 206319

#extract regionB
cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"


    #cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 459517 && $9 >= 349846 && $10 > 1000' | awk '$8 <= 349900 || $9 >= 459450' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk8 | head -n1

    startData=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 459517 && $9 >= 349846 && $10 > 1000 && $8 <= 349900' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nrk10 | head -n1 | awk '{print $3,$5}')

    regionBStart=$(echo $startData | awk '{print $1}')
    echo $regionBStart

    orientation1=$(echo $startData | awk '{print $2}')

    endData=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 459517 && $9 >= 349846 && $10 > 1000 && $9 >= 459450' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nrk10 | head -n1 | awk '{print $4,$5}')

    regionBEnd=$(echo $endData | awk '{print $1}')
    echo $regionBEnd

    orientation2=$(echo $endData | awk '{print $2}')


    if [[ $orientation1 == '-' || $orientation2 == '-' ]]; then
        echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionB/reverseOrientedMatches.txt
    fi

    #write region B coords to bed file
    #echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName}
    echo -e $chrName"\t"$regionBStart"\t"$regionBEnd"\t"regionB_${hapName} > /data/wrayva/output/extract_regions/regionB/bed/regionB_${file%.fna}.bed

    #extract with cut
    cutStart="$(($regionBStart + 1))"
    zeroIndexedStart="$(($regionBStart))"
    cutEnd="$(($regionBEnd))"
    zeroIndexedEnd=$cutEnd
    echo ">${file%.fna}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/regionB/sequences/regB_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/regionB/sequences/regB_${file}

    #output region B length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/regionB/lengths.txt
done

#reverse orientation matching region B
distal_HG04184_chr14_haplotype2-0000060
distal_HG03831_chr14_haplotype1-0000003
distal_HG02486_chr13_haplotype2-0000089

cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | grep HG00642	| grep chr13 | grep haplotype2-0000067 | awk '$8 <= 459517 && $9 >= 349846 && $10 > 1000 && $9 >= 459450' | cut -c-100 | sort -nrk10

cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | grep HG00642	| grep chr13 | grep haplotype2-0000067 | awk '$8 <= 459517 && $9 >= 349846 && $10 > 1000 && $8 <= 349900' | cut -c-100 | sort -nrk10

#check region B
#my sequence
cd /data/wrayva/output/extract_regions/regionB/
cat sequences/regB_distal_HG00099_chr14_haplotype1-0000022.fna | cut -c-100
>distal_HG03453_chr13_haplotype2-0000171:329134-438809
CTCTGCTGAGTACAGCAGGATAAAGTCCTTCATGACCCTTGGACTTTTTTATTGGGATTATTAAAAATCAGATTTCAGTATAAAAAACACCATAATTGAT
>distal_HG00099_chr13_haplotype1-0000019:269500-378798
GGGAGCATTTGGTCACTGCTCTGCTGAGTACAGCAGGATAAAGTCCTTCATGACCCTTGGACTTTTTTATTGGGATTATTAAAAATCAGATTTCAGTATA
#end
cat sequences/regB_distal_HG03453_chr13_haplotype2-0000171.fna | rev | cut -c-100
CTTTAGGAGTTATTTTATGATGGTTTGGCTTAGGTCGTCGTTTAGTTTTTCGAATAGGTGGTACTAGTTTACCCGACGTAGGGACCCTACGTTCCGACCA
#chm13 chr22 349847-459516
#beginning
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | cut -c349848-349947
GGGAGCATTTGGTCACTGCTCTGCTGAGTACAGCAGGATAAAGTCCTTCATGACCCTTGGACTTTTTTATTGGGATTATTAAAAATCAGATTTCAGTATA
#chm13 chr22
#end
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | cut -c459417-459516 | rev
CTTTAGGAGTTATTTTATGATGGTTTGGCTTAGGTCGTCGTTTAGTTTTTCGAATAGGTGGTACTAGTTTACCCGACGTAGGGACCCTACGTTCCGACCA

#check how many times the 10 char string occurs
grep -o 'GGGAGCATTT' /data/wrayva/output/sequences/distal_HG03453_chr13_haplotype2-0000171.fna | wc -l

cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo $file
    #regionB
    #AGGCAGGCTGACAGC
    grep -o 'AGGCAGGCTGA' $file | wc -l
    #GGATGTTACACCTGC
    grep -o 'GGATGTTACAC' $file | wc -l

    #regionB
    #grep -o 'GGGAGCATTT' $file | wc -l
    #grep -o 'CTTTAGGAGT' $file | wc -l
    #DJ
    #AATGTAAGGCTGAAA
    grep -o 'AATGTAAGGCT' $file | wc -l
    #CGGAAGTTCAGGTCG
    grep -o 'CGGAAGTTCAG' $file | wc -l
done


#good and bad regionB
for file in `ls /data/wrayva/output/extract_regions/regionB/sequences/regB_distal_*fna`; do
    if [[ $(cat $file | head -n2 | tail -n1 | cut -c-10) == 'GGGAGCATTT' || $(cat $file | head -n2 | tail -n1 | rev | cut -c-10) == 'CTTTAGGAGT' ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionB/good_regionB.txt
    else
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionB/bad_regionB.txt
    fi
done

for file in `ls /data/wrayva/output/extract_regions/regionB/sequences/regB_distal_*fna`; do
    file2=${file#/data/wrayva/output/extract_regions/regionB/sequences/regB_}
    echo $file2
    genomeName=$(echo ${file2%.fna} | cut -d '_' -f2)
    echo "genomeName: $genomeName"
    chrName=$(echo ${file2%.fna} | cut -d '_' -f3)
    echo "chrName: $chrName"
    hapName=$(echo ${file2%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"

    length=$(cat /data/wrayva/output/extract_regions/regionB/lengths.txt | grep $genomeName | grep $chrName | grep $hapName | awk '{print $4}')
    echo $length

    if [[ $length -ge 109000 &&  $length -le 169000 ]]; then
        echo "good"
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionB/good_regionB.txt
    else
        echo "bad"
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionB/bad_regionB.txt
    fi
done

file=/data/wrayva/output/extract_regions/regionB/sequences/regB_distal_HG04199_chr14_haplotype1-0000016.fna

module load seqtk
for name in `cat /data/wrayva/output/extract_regions/regionB/reverseOrientedMatches.txt`; do
    mv /data/wrayva/output/extract_regions/regionB/sequences/regB_${name}.fna /data/wrayva/output/extract_regions/regionB/sequences/preReverse_regB_${name}.fna
    seqtk seq -r /data/wrayva/output/extract_regions/regionB/sequences/preReverse_regB_${name}.fna  > /data/wrayva/output/extract_regions/regionB/sequences/regB_${name}.fna
done

cat preReverse_regB_distal_HG02486_chr13_haplotype2-0000089.fna | cut -c-100
cat regB_distal_HG02486_chr13_haplotype2-0000089.fna | rev | cut -c-100

Use this to reverse complement a fasta file:
seqtk seq -r in.fq > out.fq

#put all good into a fasta
for name in `cat /data/wrayva/output/extract_regions/regionB/good_regionB.txt`; do
    cat /data/wrayva/output/extract_regions/regionB/sequences/regB_${name}.fna >> /data/wrayva/output/extract_regions/regionB/good_regionB.fna
done


#run mafft
sbatch --time=30:00:00 --mem=32g --cpus-per-task=20 /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/regionB/good_regionB.fna /data/wrayva/output/extract_regions/regionB/good_regionB_alignment_mafft.aln
#30504025
cat /data/wrayva/output/slurm-30504025.out
#30790798
cat /home/wrayva/slurm-30520381.out
cat /data/wrayva/output/extract_regions/regionB/slurm-30943627.out

#run iqtree
sbatch --time=8:00:00 --mem=16g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/regionB/good_regionB_alignment_mafft.aln /data/wrayva/output/extract_regions/regionB/iqtree_output/
#30949413
cat /data/wrayva/output/extract_regions/regionB/slurm-30949413.out

cp /data/wrayva/output/extract_regions/regionB/good_regionB_alignment_mafft.aln.treefile /data/wrayva/output/extract_regions/regionB/good_regionB_alignment_mafft_copy.aln.treefile

#create new plot with super populations
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/regionB/good_regionB_alignment_mafft_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done


module load R
R

library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionB/good_regionB_alignment_mafft_copy.aln.treefile")
chrs <- list(chr13=c(tree$tip.label[grep("chr13", tree$tip.label)]),
             chr14=c(tree$tip.label[grep("chr14", tree$tip.label)]),
             chr15=c(tree$tip.label[grep("chr15", tree$tip.label)]),
             chr21=c(tree$tip.label[grep("chr21", tree$tip.label)]),
             chr22=c(tree$tip.label[grep("chr22", tree$tip.label)])
)
grouped_tree1 <- groupOTU(tree, chrs,
                       group_name = "chromosome")
ggtree(grouped_tree1, aes(color=chromosome), layout='circular') +
  theme(legend.position=c(.95,.95), legend.key.size = unit(.1, 'cm')) +
  geom_tiplab(align=TRUE, linesize=.1, size=1)

ggsave("/data/wrayva/output/plots/regionBTree_circular_color_by_chr.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

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
ggsave("/data/wrayva/output/plots/regionBTree_circular_color_by_super_pop.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)



#get pop statistics
cd /data/wrayva/output/extract_regions/regionB
eur_count=0
eas_count=0
amr_count=0
afr_count=0
sas_count=0
#for file in `ls distal_*fna`; do
for file in `cat good_regionB.txt`; do
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
#regionB: EUR:23 EAS:58 AMR:86 AFR:91 SAS:68


#compare chr21 and chr22 to make sure we extract correct regB
${scriptDir}/run_minimap.sh /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.fa /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr21.distal.fa $regionDir/chm13chr22vschr21.paf

${scriptDir}/run_minimap.sh $regionDir/chm13chr22regionB.fa /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr21.distal.fa $regionDir/chm13chr22regBvschr21.paf

#try making region B the Query
${scriptDir}/run_minimap.sh /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr21.distal.fa $regionDir/chm13chr22regionB.fa $regionDir/chm13chr22regBqueryvschr21ref.paf

#filter the paf
awk '($17~"tp:A:P")  {print $0}' chm13chr22vschr21.paf > filtered_chm13chr22vschr21.paf

#chr22:1-4793795
chr21:1-3132131
sed -i "s/chr21:1-3132131/chr21/g" $regionDir/chm13chr22regBqueryvschr21ref.paf
sed -i "s/chr22:226744-459516/chr22/g" $regionDir/chm13chr22regBqueryvschr21ref.paf
