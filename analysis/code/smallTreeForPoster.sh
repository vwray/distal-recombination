#pick 10 random haps from each superpop
#cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
cd /data/wrayva/output/extract_regions/regionA2
#for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
for file in `cat good_regionA.txt`; do
    echo ${file}
    genomeName=$(echo ${file} | cut -d '_' -f2)
    chrName=$(echo ${file} | cut -d '_' -f3)
    hapName=$(echo ${file} | cut -d '_' -f4)
    echo "hapName: $hapName"
    #echo $genomeName
    #cat /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_copy2.fna.treefile | grep -o distal_$genomeName | wc -l

    #file=/data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2_COPY.aln.treefile
    #toReplace=distal_$genomeName
    #echo $toReplace

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    echo $file >> /data/wrayva/output/extract_regions/small_regionA/${superpop}_haps.txt

    #replaceWith=${superpop}_${pop}_${genomeName}
    #echo $replaceWith

    #sed -i "s/$toReplace/$replaceWith/g" $file

    #write to files for each superpop
done

cd /data/wrayva/output/extract_regions/small_regionA
for hapFile in `ls *_haps.txt`; do
    for file in `sort -R $hapFile | head -n10`; do
        echo $file >> ten_haps_from_each_pop.txt
    done
done




cd /data/wrayva/output/extract_regions/small_regionA
for file in `cat ten_haps_from_each_pop.txt`; do
    echo ${file}
    genomeName=$(echo ${file} | cut -d '_' -f2)
    chrName=$(echo ${file} | cut -d '_' -f3)
    hapName=$(echo ${file} | cut -d '_' -f4)
    echo "hapName: $hapName"
done


#put all random samples in one fasta
cd /data/wrayva/output/extract_regions/small_regionA
for name in `cat ten_haps_from_each_pop.txt`; do
    cat /data/wrayva/output/extract_regions/regionA2/sequences/regA_${name}.fna >> /data/wrayva/output/extract_regions/small_regionA/random_samples.fna
done




#run mafft
sbatch --time=8:00:00 --mem=16g --cpus-per-task=10 /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/small_regionA/random_samples.fna /data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft.aln
#31004726
cat /data/wrayva/output/extract_regions/small_regionA/slurm-31004726.out

#run iqtree
sbatch --time=8:00:00 --mem=16g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft.aln /data/wrayva/output/extract_regions/small_regionA/iqtree_output/
#30949413
cat /data/wrayva/output/extract_regions/small_regionA/slurm-30949413.out




#create new plot with super populations
cp /data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft.aln.treefile /data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft_copy.aln.treefile

cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done

cd /data/wrayva/output/extract_regions/small_regionA

file=/data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft_copy.aln.treefile
for hapString in `grep -o 'haplotype[1-2]-[0-9]\+_[0-9]\+-[0-9]\+' $file`; do
    toReplace=$hapString
    echo $toReplace
    replaceWith=hap$(echo $hapString | cut -c 10)
    echo $replaceWith
    sed -i "s/$toReplace/$replaceWith/g" $file
done




theme_tree2(plot.margin=margin(10, 10, 10, 10)) +





module load R
R

library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft_copy.aln.treefile")
chrs <- list(chr13=c(tree$tip.label[grep("chr13", tree$tip.label)]),
             chr14=c(tree$tip.label[grep("chr14", tree$tip.label)]),
             chr15=c(tree$tip.label[grep("chr15", tree$tip.label)]),
             chr21=c(tree$tip.label[grep("chr21", tree$tip.label)]),
             chr22=c(tree$tip.label[grep("chr22", tree$tip.label)])
)
grouped_tree1 <- groupOTU(tree, chrs,
                       group_name = "chromosome")
pdf(file = "/data/wrayva/output/plots/small_regionATree_circular_color_by_chr2.pdf")
ggtree(grouped_tree1, aes(color=chromosome), size=3, layout='circular') +
  theme(legend.position=c(1.0,.98), legend.key.size = unit(2, 'cm'),plot.margin=margin(150, 150, 150, 150), legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(2))) +
  geom_tiplab(align=TRUE, linesize=.5, size=6, fontface='bold')

ggsave("/data/wrayva/output/plots/small_regionATree_circular_color_by_chr_widthheight=60size=6boldmargin150theme2legendfontsize2.pdf", width = 60, height = 60, units = "cm", limitsize = FALSE)

pops <- list(EUR=c(tree$tip.label[grep("EUR", tree$tip.label)]),
             EAS=c(tree$tip.label[grep("EAS", tree$tip.label)]),
             AMR=c(tree$tip.label[grep("AMR", tree$tip.label)]),
             AFR=c(tree$tip.label[grep("AFR", tree$tip.label)]),
             SAS=c(tree$tip.label[grep("SAS", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, pops,
                       group_name = "population")
ggtree(grouped_tree2, aes(color=population), size=3, layout='circular') +
 theme(legend.position=c(1.0,.98), legend.key.size = unit(2, 'cm'),plot.margin=margin(150, 150, 150, 150), legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(2))) +
 geom_tiplab(align=TRUE, linesize=.5, size=6, fontface='bold')

ggsave("/data/wrayva/output/plots/small_regionATree_circular_color_by_pop_widthheight=60size=6boldmargin150theme2legendfontsize2.pdf", width = 60, height = 60, units = "cm", limitsize = FALSE)











cd /data/wrayva/output/extract_regions/small_regionA
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy.aln.mldist
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done

cd /data/wrayva/output/extract_regions/small_regionA
for line in `cat /data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2_copy.aln.mldist | tail n$(cat small_regionA_alignment_mafft_copy.aln.mldist | head -n1)`; do
    awk -F '[[:blank:],]' '{m=$1; for (i=2; i<=NF; i++) if ($i < m) m = $i; print m}'
done


cat small_regionA_alignment_mafft_copy.aln.mldist | head -n2 | tail -n1 | awk -F '[[:blank:],]' '{m=$1; for (i=2; i<=NF; i++) if ($i < m) m = $i; print m}'



python


import numpy as np
import random
import pandas as p
#file_obj = open("/data/wrayva/output/extract_regions/DJ/good_dj_alignment_mafft_copy.aln.mldist", "r")
file_obj = open("/data/wrayva/output/extract_regions/small_regionA/small_regionA_alignment_mafft_copy.aln.mldist", "r")

file_data = file_obj.read()
lines = file_data.splitlines()
matrix=[[0 for i in range(len(lines))] for j in range(len(lines))]
for i in range(0,len(lines)):
    print(lines[i].split())
    matrix[i]=lines[i].split()


#find min in each row, ignoring first col and the diag [i+1][i]

#for row in matrix:
countOfPopMatch=0
countOfChrMatch=0
for i in range(len(matrix)):
    if(i==0):
        continue
    row=matrix[i]
    print(row)
    #print(row[1:])
    #print(p.Series(row[1:]))
    #index = p.Series(row[1:]).idxmin()
    minimum=min(row[1:])
    #print(index, "is the min")
    currentMin=2.0
    currentMinIndex=-1
    for j in range(len(row)):
        if(j!=0 and j!=i and float(row[j])!=0.0 and float(row[j])<currentMin):
            currentMin=float(row[j])
            currentMinIndex=j
    print(f'The min for row {i} is {currentMin} at column {currentMinIndex}')
    #check if pop matches
    thisPop=row[0][0:3]
    neighborPop=matrix[j][0][0:3]
    print(f'pop of this row is {row[0][0:3]}')
    print(f'pop of closest neighbor is {matrix[j][0][0:3]}')
    if(thisPop==neighborPop):
        print("same pop")
        countOfPopMatch+=1
    thisChr=row[0][16:21]
    neighborChr=matrix[j][0][16:21]
    print(thisChr)
    print(neighborChr)
    if(thisChr==neighborChr):
        print("same chr")
        countOfChrMatch+=1

#for small tree, 10/50=.192 for same pop, 7/50=0.14 for same chr
#for region A2 tree, 60/272=.2206 for same pop, 53/272=.1949 for same chr
#for DJ tree,65/267=0.2434 have their nearest neighbor as same pop and 50/267=0.1873 have their nearest neighbor as same chr


print(f'for small region A tree,{countOfPopMatch}/{len(matrix)-1}={float(countOfPopMatch)/float(len(matrix)-1)} have their nearest neighbor as same pop and {countOfChrMatch}/{len(matrix)-1}={float(countOfChrMatch)/float(len(matrix)-1)} have their nearest neighbor as same chr')










x=[]
for i in range(0, len(lines)):
    x.append(int(lines[i]))
