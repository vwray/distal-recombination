#extract region C
cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"


    #cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 459517 && $9 >= 349846 && $10 > 1000' | awk '$8 <= 349900 || $9 >= 459450' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk8 | head -n1

    startData=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 3841140 && $9 >= 3634819 && $10 > 1000 && $8 <= 3634900 && $8 >= 3634750' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nrk10 | head -n1 | awk '{print $3,$5}')

    regionCStart=$(echo $startData | awk '{print $1}')
    echo $regionCStart

    orientation1=$(echo $startData | awk '{print $2}')

    endData=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 3841140 && $9 >= 3634819 && $10 > 1000 && $9 >= 3841100 && $9 <= 3841200' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nrk10 | head -n1 | awk '{print $4,$5}')

    regionCEnd=$(echo $endData | awk '{print $1}')
    echo $regionCEnd

    orientation2=$(echo $endData | awk '{print $2}')


    if [[ $orientation1 == '-' || $orientation2 == '-' ]]; then
        echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionC/reverseOrientedMatches.txt
    fi

    #write region C coords to bed file
    #echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName}
    echo -e $chrName"\t"$regionCStart"\t"$regionCEnd"\t"regionC_${hapName} > /data/wrayva/output/extract_regions/regionC/bed/regionC_${file%.fna}.bed

    #extract with cut
    cutStart="$(($regionCStart + 1))"
    zeroIndexedStart="$(($regionCStart))"
    cutEnd="$(($regionCEnd))"
    zeroIndexedEnd=$cutEnd
    echo ">${file%.fna}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/regionC/sequences/regC_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/regionC/sequences/regC_${file}

    #output region C length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/regionC/lengths.txt
done





for file in `ls /data/wrayva/output/extract_regions/regionC/sequences/regC_distal_*fna`; do
    file2=${file#/data/wrayva/output/extract_regions/regionC/sequences/regC_}
    echo $file2
    genomeName=$(echo ${file2%.fna} | cut -d '_' -f2)
    echo "genomeName: $genomeName"
    chrName=$(echo ${file2%.fna} | cut -d '_' -f3)
    echo "chrName: $chrName"
    hapName=$(echo ${file2%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"

    length=$(cat /data/wrayva/output/extract_regions/regionC/lengths.txt | grep $genomeName | grep $chrName | grep $hapName | awk '{print $4}')
    echo $length

    if [[ $length -ge 205000 &&  $length -le 207000 ]]; then
        echo "good"
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionC/good_regionC.txt
    else
        echo "bad"
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionC/bad_regionC.txt
    fi
done

file=/data/wrayva/output/extract_regions/regionC/sequences/regC_distal_HG04199_chr14_haplotype1-0000016.fna

module load seqtk
for name in `cat /data/wrayva/output/extract_regions/regionC/reverseOrientedMatches.txt`; do
    mv /data/wrayva/output/extract_regions/regionC/sequences/regC_${name}.fna /data/wrayva/output/extract_regions/regionC/sequences/preReverse_regC_${name}.fna
    seqtk seq -r /data/wrayva/output/extract_regions/regionC/sequences/preReverse_regC_${name}.fna  > /data/wrayva/output/extract_regions/regionC/sequences/regC_${name}.fna
done

cat preReverse_regC_distal_HG02486_chr13_haplotype2-0000089.fna | cut -c-100
cat regC_distal_HG02486_chr13_haplotype2-0000089.fna | rev | cut -c-100

Use this to reverse complement a fasta file:
seqtk seq -r in.fq > out.fq

#put all good into a fasta
for name in `cat /data/wrayva/output/extract_regions/regionC/good_regionC.txt`; do
    cat /data/wrayva/output/extract_regions/regionC/sequences/regC_${name}.fna >> /data/wrayva/output/extract_regions/regionC/good_regionC.fna
done


#run mafft
sbatch --time=30:00:00 --mem=32g --cpus-per-task=20 /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/regionC/good_regionC.fna /data/wrayva/output/extract_regions/regionC/good_regionC_alignment_mafft.aln
#30945515
cat /data/wrayva/output/extract_regions/regionC/slurm-30945515.out

#run iqtree
sbatch --time=8:00:00 --mem=16g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/regionC/good_regionC_alignment_mafft.aln /data/wrayva/output/extract_regions/regionC/iqtree_output/
#30949538
cat /data/wrayva/output/extract_regions/regionC/slurm-30949538.out



cp /data/wrayva/output/extract_regions/regionC/good_regionC_alignment_mafft.aln.treefile /data/wrayva/output/extract_regions/regionC/good_regionC_alignment_mafft_copy.aln.treefile

#create new plot with super populations
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/regionC/good_regionC_alignment_mafft_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done

cd /data/wrayva/output/extract_regions/regionC/


module load R
R

library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionC/good_regionC_alignment_mafft_copy.aln.treefile")
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

ggsave("/data/wrayva/output/plots/regionCTree_circular_color_by_chr.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

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
ggsave("/data/wrayva/output/plots/regionCTree_circular_color_by_super_pop.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)



#get pop statistics
cd /data/wrayva/output/extract_regions/regionC
eur_count=0
eas_count=0
amr_count=0
afr_count=0
sas_count=0
#for file in `ls distal_*fna`; do
for file in `cat good_regionC.txt`; do
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
#region C: EUR:20 EAS:46 AMR:62 AFR:47 SAS:58
