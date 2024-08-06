cd /data/wrayva/output/extract_regions
#bam file
samtools view -L regionA.bed /data/wrayva/output/minimap_chm13chr22masked/filtered.bam | wc -l

samtools view -L regionA.bed /data/wrayva/output/minimap_chm13chr22masked/filtered.bam | cut -c-200 | grep distal_HG00099_chr13_haplotype1-0000019 | sort -nk8

#paf file
cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_aln.paf | awk '$8 <= 67033 && $9 >= 4577' | grep distal_HG00099_chr13_haplotype1-0000019 | sort -k1 | head -n100

cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577' | grep distal_HG00099_chr13_haplotype1-0000019



#pick largest
#double check that it's reasonable
#possibly pick top 2 largest, if one is invert

#or, take smallest start and largest end?

#regionA
cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"


    #cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk8 | head -n1


    regionAStart=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577 && $3 < 100000 && $10 > 1000 && $8 <= 4600' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk3 | head -n1 | awk '{print $3}')
    echo $regionAStart

    regionAEnd=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577 && $3 < 100000 && $10 > 1000 && $9 >= 67000' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk4 | head -n1 | awk '{print $4}')
    echo $regionAEnd


    if [[ $orientation1 == '-' || $orientation2 == '-' ]]; then
        echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionA2/reverseOrientedMatches.txt
    fi

    #write region A coords to bed file
    #echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName}
    echo -e $chrName"\t"$regionAStart"\t"$regionAEnd"\t"regionA_${hapName} > /data/wrayva/output/extract_regions/regionA2/regionA_${file%.fna}.bed

    #extract with cut
    cutStart="$(($regionAStart + 1))"
    zeroIndexedStart="$(($regionAStart))"
    cutEnd="$(($regionAEnd))"
    zeroIndexedEnd=$cutEnd
    echo ">${file%.fna}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/regionA2/sequences/regA_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/regionA2/sequences/regA_${file}

    #output region A length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/regionA2/lengths.txt
done

cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577 && $3 < 100000 && $10 > 1000' | grep distal_HG02280_chr21_haplotype1-0000018 | sort -nk3



#check region A
#my sequence
cd /data/wrayva/output/extract_regions/regionA/
cat sequences/regA_distal_HG03453_chr13_haplotype2-0000171.fna | cut -c-100
AGGCAGGCTGACAGCGGTCATGTTTCTGCCTACAGCGCCTGCCTATCTCTTTTGAATGTCCTTCTCTACCCTACTCTGTACTTATGGTGCCAGGTTTCTC
#end
cat sequences/regA_distal_HG03453_chr13_haplotype2-0000171.fna | rev | cut -c-100
GGATGTTACACCTGCTTAGAACTTTTATAACACGATGCACAATACAGTCAGTATTGTCGAGGGGATAACACGTTAAGGCGAATATCCTTTACAGGTGTTA
#chm13 chr22
#beginning
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | cut -c4579-4678
AGGCAGGCTGACAGCGGTCAAGTTTCTGCCTACAGCGCCTGCCTATCTCTTTTGAATGTCCTTCTCTACCCTACTCTGTACTTATGGTGCCAGGTTTCTC
#chm13 chr22
#end
cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa | cut -c66933-67032 | rev
GGATGTTACACCTGCTTAGAACTTTTATAACACGATGCACAATACAGTCAGTATTGTCGAGGGGATAACACGTTAAGGCGAATATCCTTTACAGGTGTTA

#ones that look different:
>distal_NA20805_chr13_haplotype2-0000136:3277-53128
AACCTACAGTTCACTGTGGTGAAGGAGCAATTCAGACTCACCCTGCAGCTTCCTGTATAGGTCAAACAGTGGCCGCTCCCTAAGGAAAAGTTGGGGGTGC
>distal_NA20752_chr15_haplotype1-0000017:17745-80769
TCACTCCCGGCCAGAGAAAAACCCTGTAGACATTAGTCACTCCTCATTCTGTCTCAAACACTCTCCCTGACCCTCAGCCCTAGGTAGCAACTACCTAGTG
>distal_NA20752_chr13_haplotype1-0000003:11877-100073
TCACTCCCGGCCAGAGAAAAACCCTGTAGACATTAGTCACTCCTCATTCTGTCTCAAACACTCTCCCTGACCCTCAGCCCTAGGTAGCAACTACCTAGTG
>distal_NA18983_chr21_haplotype2-0000056:5150-64701
AACCAACCCAAATATCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATACACCATGGAATACTATGCAGCCATAAAAAATGAAGAGTTCATGT
(more in /data/wrayva/output/extract_regions/regionA/sequences/abnormal_regionA.txt)


#check on sequences with region A missing



#bad regionA
for file in `ls /data/wrayva/output/extract_regions/regionA2/sequences/regA_distal_*fna`; do
    if [[ $(cat $file | head -n2 | tail -n1 | cut -c-10) != 'AGGCAGGCTG' ]]; then
        cut -c-100 $file
    fi
done

#good and bad regionA
for file in `ls /data/wrayva/output/extract_regions/regionA2/sequences/regA_distal_*fna`; do
    if [[ $(cat $file | head -n2 | tail -n1 | cut -c-10) == 'AGGCAGGCTG' && $(cat $file | head -n2 | tail -n1 | rev | cut -c-10) == 'GGATGTTACA' ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionA2/good_regionA.txt
    else
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionA2/bad_regionA.txt
    fi
done

#put all good ones in one fasta
for name in `cat /data/wrayva/output/extract_regions/regionA2/good_regionA.txt`; do
    cat /data/wrayva/output/extract_regions/regionA2/sequences/regA_${name}.fna >> /data/wrayva/output/extract_regions/regionA2/good_regionA.fna
done

# put all good ones in separate fastas
for name in `cat /data/wrayva/output/extract_regions/regionA2/good_regionA.txt`; do
    cp /data/wrayva/output/extract_regions/regionA2/sequences/regA_${name}.fna /data/wrayva/output/extract_regions/regionA2/good_sequences
done

#run Muscle
sbatch --time=30:00:00 --mem=256g --cpus-per-task=20 /data/wrayva/scripts/run_muscle.sh align /data/wrayva/output/extract_regions/regionA/good_regionA.fna /data/wrayva/output/extract_regions/regionA/good_regionA_muscle_alignment.fna
#30494115

#run mafft
sbatch --time=30:00:00 --mem=256g --cpus-per-task=20 /data/wrayva/scripts/run_mafft.sh /data/wrayva/output/extract_regions/regionA/good_regionA.fna /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft.aln
#30500980
cat ./extract_regions/regionA/slurm-30500980.out
#It worked!

#run mafft with retree 2
sbatch --time=30:00:00 --mem=128g --cpus-per-task=20 /data/wrayva/scripts/run_mafft_retree2.sh /data/wrayva/output/extract_regions/regionA2/good_regionA.fna /data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2.aln
#30759224

#run parsnp
sbatch --time=10:00:00 --mem=32g /data/wrayva/scripts/run_parsnp.sh /data/wrayva/output/extract_regions/regionA2/chr22_regionA.fna /data/wrayva/output/extract_regions/regionA2/good_sequences /data/wrayva/output/extract_regions/regionA2/parsnp-output -c

#run iqtree (create script)
touch /data/wrayva/scripts/run_iqtree.sh

sbatch --time=8:00:00 --mem=32g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2.aln /data/wrayva/output/extract_regions/regionA2/iqtree_output/
#30765070
cat /data/wrayva/output/extract_regions/regionA/slurm-30508028.out

output from iqtree: /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft.fna.treefile


#view in R:
module load R
R

# It works! Colors by chromosome.
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2_COPY.aln.treefile")
pdf(file = "/data/wrayva/output/plots/regionATree_circular_color_retree2_1.pdf")
chrs <- list(chr13=c(tree$tip.label[grep("chr13", tree$tip.label)]),
             chr14=c(tree$tip.label[grep("chr14", tree$tip.label)]),
             chr15=c(tree$tip.label[grep("chr15", tree$tip.label)]),
             chr21=c(tree$tip.label[grep("chr21", tree$tip.label)]),
             chr22=c(tree$tip.label[grep("chr22", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, chrs,
                       group_name = "chromosome")
ggtree(grouped_tree2, aes(color=chromosome), layout='circular') +
  theme(legend.position=c(.95,.95), legend.key.size = unit(.1, 'cm')) +
  geom_tiplab(align=TRUE, linesize=.1, size=1)
dev.off()

ggsave("/data/wrayva/output/plots/regionA2Tree_circular_color_by_chr.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

#, branch.length='none' goes inside ggtree

#color by population
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_retree2_COPY.aln.treefile")
pdf(file = "/data/wrayva/output/plots/regionATree_circular_color_by_pop_retree2_1.pdf")
pops <- list(GBR=c(tree$tip.label[grep("GBR", tree$tip.label)]),
             FIN=c(tree$tip.label[grep("FIN", tree$tip.label)]),
             CHS=c(tree$tip.label[grep("CHS", tree$tip.label)]),
             PUR=c(tree$tip.label[grep("PUR", tree$tip.label)]),
             CLM=c(tree$tip.label[grep("CLM", tree$tip.label)]),
             ACB=c(tree$tip.label[grep("ACB", tree$tip.label)]),
             PEL=c(tree$tip.label[grep("PEL", tree$tip.label)]),
             KHV=c(tree$tip.label[grep("KHV", tree$tip.label)]),
             CDX=c(tree$tip.label[grep("CDX", tree$tip.label)]),
             GWD=c(tree$tip.label[grep("GWD", tree$tip.label)]),
             PJL=c(tree$tip.label[grep("PJL", tree$tip.label)]),
             ESN=c(tree$tip.label[grep("ESN", tree$tip.label)]),
             MSL=c(tree$tip.label[grep("MSL", tree$tip.label)]),
             BEB=c(tree$tip.label[grep("BEB", tree$tip.label)]),
             STU=c(tree$tip.label[grep("STU", tree$tip.label)]),
             YRI=c(tree$tip.label[grep("YRI", tree$tip.label)]),
             CHB=c(tree$tip.label[grep("CHB", tree$tip.label)]),
             JPT=c(tree$tip.label[grep("JPT", tree$tip.label)]),
             LWK=c(tree$tip.label[grep("LWK", tree$tip.label)]),
             TSI=c(tree$tip.label[grep("TSI", tree$tip.label)]),
             GIH=c(tree$tip.label[grep("GIH", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, pops,
                       group_name = "population")
ggtree(grouped_tree2, aes(color=population), layout='circular') +
  theme(legend.position=c(.9,.9), legend.key.size = unit(.2, 'cm')) +
  geom_tiplab(align=TRUE, linesize=.1, size=1)
dev.off()

#color by super population
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2_COPY.aln.treefile")
pdf(file = "/data/wrayva/output/plots/regionA2Tree_circular_color_by_super_pop_retree2_1.pdf")
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
  geom_tiplab(align=TRUE, linesize=.1, size=3)
dev.off()


ggsave("/data/wrayva/output/plots/regionA2Tree_circular_color_by_super_pop_retree2_3.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

#+coord_cartesian(clip="off")
#ggplot2::xlim(0, .4)


#try msa plot to visualize msa, look for gaps
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_retree2.aln.treefile")
pdf(file = "/data/wrayva/output/plots/regionATree_cmsa_retree2_2.pdf")

msaplot(p=ggtree(tree), fasta="/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_retree2_COPY.aln")
dev.off()



msaplot(p=ggtree(tree), fasta="/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_retree2_COPY.aln", window=c(150, 175))


/data/wrayva/scripts/plot_msa.r

#!/bin/bash
module load R
R --no-echo --no-restore --no-save < /data/wrayva/scripts/plot_msa.r > /data/wrayva/output/plots/msa_r_job.out

sbatch --time=8:00:00 --mem=32g plot_msa_with_r.sh
#30705456

#get pop label for each genome, replace in tree file
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    echo $genomeName
    #cat /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_COPY.fna.treefile | grep -o distal_$genomeName | wc -l

    file=/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_retree2_COPY.aln.treefile
    toReplace=distal_$genomeName
    echo $toReplace

    replaceWith=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)_${genomeName}
    echo $replaceWith

    sed -i "s/$toReplace/$replaceWith/g" $file
done

sed -i "s/:/_/g" /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_retree2_COPY.aln


#create new plot with super populations
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    #echo $genomeName
    #cat /data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_copy2.fna.treefile | grep -o distal_$genomeName | wc -l

    file=/data/wrayva/output/extract_regions/regionA2/good_regionA_alignment_mafft_retree2_COPY.aln.treefile
    toReplace=distal_$genomeName
    #echo $toReplace

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}
    #echo $replaceWith

    sed -i "s/$toReplace/$replaceWith/g" $file
done


#create moddotplots
source myconda
conda activate dfv
cd /data/wrayva/gitRepos/ModDotPlot
moddotplot static -f /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa /data/wrayva/output/sequences/distal_HG02280_chr21_haplotype1-0000018.fna -o /data/wrayva/output/moddotplot/gap_in_alignment_HG02280_chr21_haplotype1-0000018 --compare-only

#extract chr22 region A
cd /data/wrayva/output/extract_regions/regionA/moddotplot_test/
#echo ">chr22:1-67032" > chr22_regionA.fna
#cat /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa
module load samtools
samtools faidx /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa chr22:1-67032 > chr22_regionA.fna

samtools faidx /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa chr22:4578-67032 > chr22_regionA.fna

samtools faidx /data/wrayva/output/sequences/distal_HG02280_chr21_haplotype1-0000018.fna distal_HG02280_chr21_haplotype1-0000018:1-100000 > distal_HG02280_chr21_haplotype1-0000018_first100000.fna



cd /data/wrayva/gitRepos/ModDotPlot
moddotplot static -f /data/wrayva/output/extract_regions/regionA/moddotplot_test/chr22_regionA.fna /data/wrayva/output/extract_regions/regionA/moddotplot_test/distal_HG02280_chr21_haplotype1-0000018_first100000.fna -o /data/wrayva/output/moddotplot/regionA_gap_in_alignment_HG02280_chr21_haplotype1-0000018 --compare-only

#trim MSA
samtools faidx good_regionA_alignment_mafft.aln *:1-10 > trim_msa/chr22_regionA.fna



#chm13 chr4 extract region 5955-68255
#chm13 chr4 gzipped: /data/pickettbd/parent-of-origin/chm13/v2.0/split-per-chr/chr4.fa.gz
#780*80=62,400
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4
module load samtools
samtools faidx /data/pickettbd/parent-of-origin/chm13/v2.0/split-per-chr/chr4.fa.gz chr4:5955-68255 > /data/wrayva/output/extract_regions/regionA2/chm13chr4/chm13chr4regionA.fa

#extracted region A here:
/data/wrayva/output/extract_regions/regionA2/chm13chr4/chm13chr4regionA.fa

#minimap map this to all Assemblies, see if present on chromosomes outside of acros and if found on chr4 consistently
module load minimap2
loop through assembly files, for each one, run minimap with chm13chr4regionA as ref

minimap2 -x asm20 /data/wrayva/output/sequences/distal_NA20905_chr22_haplotype2-0000055.fna /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa > /data/wrayva/output/minimap_dj/NA20905_chr22_haplotype2-0000055_ref_to_DJ_query.paf

#skip the above, for now just look at existing assembly-ref.norm.mashmap files
#first argument is the directory where genomes are located
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
#second argument is the output directory
outputDir=/data/wrayva/output/extract_regions/regionA2/chm13chr4
cd $genomeDir
#loop through genomes
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    echo $genomeName
    cd $genomeDir/$genomeName/verkko-hi-c
    #get haplotypes that contain match with the our region A of chr 4 5955-68255
    cat analysis/assembly-ref.norm.mashmap | awk '$6 == "chr4" && $8 <= 68256 && $9 >= 5954 && $8 <=  10000' | awk '{if ($(NF-1) > 0.99 && $4-$3 > 1000) print $0}' | while read -r mashmapLine; do
    #|sed s/id:f://g |awk '{if ($(NF-1) > 0.99 && $4-$3 > 200000) print $0}' | sort -nk3 | while read -r mashmapLine; do
        #look up chromosome label
        hapName=$(echo $mashmapLine | awk '{print $1}')
        chrAndHap=$(grep $hapName analysis/assembly.refOriented_reAssigned.fasta.bz.fai | grep 'chr*' | awk '{print $1}')
        if [[ ! -z "$chrAndHap" ]]; then
            echo "${genomeName}_${chrAndHap}" >> $outputDir/chm13chr4regionA_matches.txt
        fi
    done
done



genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
#second argument is the output directory
outputDir=/data/wrayva/output/extract_regions/regionA2/chm13chr4
cd $genomeDir
#loop through genomes
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    #echo $genomeName
    cd $genomeDir/$genomeName/verkko-hi-c
    #get haplotypes that contain match with the our region A of chr 4 5955-68255
    cat analysis/assembly-ref.norm.mashmap | awk '$6 == "chr4" && $8 <= 68256 && $9 >= 5954 && $8 <=  10000' | awk '{if ($(NF-1) > 0.99 && $4-$3 > 1000) print $0}' | awk '{print $11}'
    done
done


genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
outputDir=/data/wrayva/output/extract_regions/regionA2
cd $genomeDir
#loop through genomes
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    #echo $genomeName
    cd $genomeDir/$genomeName/verkko-hi-c
    #get haplotypes that contain match with the our region A of chr 4 5955-68255
    cat analysis/assembly-ref.norm.mashmap | awk '$6 == "chr4" && $8 <= 68256 && $9 >= 5954 && $8 <=  10000' | awk '{if ($(NF-1) > 0.99 && $4-$3 > 1000) print $0}' | while read -r mashmapLine; do
    #|sed s/id:f://g |awk '{if ($(NF-1) > 0.99 && $4-$3 > 200000) print $0}' | sort -nk3 | while read -r mashmapLine; do
        #look up chromosome label
        hapName=$(echo $mashmapLine | awk '{print $1}')
        chrAndHap=$(grep $hapName analysis/assembly.refOriented_reAssigned.fasta.bz.fai | grep 'chr*' | awk '{print $1}')
        if [[ ! -z "$chrAndHap" && $(echo $chrAndHap | cut -c-5) == "chr10" ]]; then
            echo "${genomeName}_${chrAndHap}"
            echo $mashmapLine >> $outputDir/chm13chr4regionA_matches.txt
        fi
    done
done

chr4matches=$(grep chr4 $outputDir/chm13chr4regionA_matches.txt | wc -l)
chr13matches=$(grep chr13 $outputDir/chm13chr4regionA_matches.txt | wc -l)
chr14matches=$(grep chr14 $outputDir/chm13chr4regionA_matches.txt | wc -l)
chr15matches=$(grep chr15 $outputDir/chm13chr4regionA_matches.txt | wc -l)
chr21matches=$(grep chr21 $outputDir/chm13chr4regionA_matches.txt | wc -l)
chr22matches=$(grep chr22 $outputDir/chm13chr4regionA_matches.txt | wc -l)
chr10matches=$(grep chr10 $outputDir/chm13chr4regionA_matches.txt | wc -l)
echo "chr4: $chr4matches chr13: $chr13matches chr14: $chr14matches chr15: $chr15matches chr21: $chr21matches chr22: $chr22matches chr10: $chr10matches"
#chr4: 201 chr13: 30 chr14: 44 chr15: 40 chr21: 42 chr22: 36 chr10: 6

$outputDir/chm13chr4regionA_matches.txt


#faidx the haplotypes in this file. create one fasta containing all these sequences
outputDir=/data/wrayva/output/extract_regions/regionA2/chm13chr4

module load samtools
cd $outputDir
for file in `cat $outputDir/chm13chr4regionA_matches.txt | sort -k1 | uniq | grep -e 'chr4' -e 'chr10'`; do
    genomeName=$(echo ${file} | cut -d '_' -f1)
    chrName=$(echo ${file} | cut -d '_' -f2)
    hapName=$(echo ${file} | cut -d '_' -f3)
    echo "${genomeName}_${chrName}_$hapName"
    #samtools faidx $genomeDir/$genomeName/verkko-hi-c/analysis/assembly.refOriented_reAssigned.fasta.bz ${chrName}_$hapName >> $outputDir/chr10and4assemblies.fa
    #grep ${chrName}_$hapName $genomeDir/$genomeName/verkko-hi-c/analysis/assembly.refOriented_reAssigned.fasta.bz.fai
    samtools faidx $genomeDir/$genomeName/verkko-hi-c/analysis/assembly.refOriented_reAssigned.fasta.bz ${chrName}_$hapName >> $outputDir/sequences/${genomeName}_${chrName}_$hapName.fa
    #toReplace=">${chrName}_$hapName"
    #echo $toReplace

    replaceWith=">${genomeName}_${chrName}_$hapName"
    echo $replaceWith

    #sed -i "s/$toReplace/$replaceWith/g" $outputDir/chr10and4assemblies.fa
    #replace first line of each file with the proper sequence identifier
    sed -i "1s/.*/$replaceWith/" $outputDir/sequences/${genomeName}_${chrName}_$hapName.fa
done

#check files to make sure they look okay
cd $outputDir/sequences
for file in `ls *.fa`; do
    cat $file | head -n2 | cut -c-100
done

#copy all fasta files to one fasta file
cd $outputDir/sequences
for file in `ls *.fa`; do
    cat $file >> $outputDir/chr10and4assemblies.fa
done

#run minimap against chm13 chr4 regA to get coords to extract
module load minimap2
minimap2 -cx asm20 /data/wrayva/output/extract_regions/regionA2/chm13chr4/chm13chr4regionA.fa /data/wrayva/output/extract_regions/regionA2/chm13chr4/chr10and4assemblies.fa > /data/wrayva/output/extract_regions/regionA2/chm13chr4/minimap_chr4_regA.paf

sbatch --time=5:00:00 --mem=32g /data/wrayva/output/extract_regions/regionA2/chm13chr4/run_minimap.sh

# filter paf file
awk '($17~"tp:A:P")  {print $0}' /data/wrayva/output/extract_regions/regionA2/chm13chr4/minimap_chr4_regA.paf > /data/wrayva/output/extract_regions/regionA2/chm13chr4/filtered_minimap_chr4_regA.paf

# extract sequences at coords in filtered paf file
# similar to below logic, but our ref is the whole thing we're interested in
module load samtools
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4/full_sequences
for file in `ls *.fa`; do
    echo ${file%.fa}
    genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
    chrName=$(echo ${file%.fa} | cut -d '_' -f2)
    hapName=$(echo ${file%.fa} | cut -d '_' -f3)
    echo "hapName: $hapName"

    coords=$(cat /data/wrayva/output/extract_regions/regionA2/chm13chr4/filtered_minimap_chr4_regA.paf | grep ${genomeName}_${chrName}_${hapName} | sort -nrk10 | head -n1 | awk '{print $3,$4,$5}')

    regionAStart=$(echo $coords | awk '{print $1}')
    echo $regionAStart

    regionAEnd=$(echo $coords | awk '{print $2}')
    echo $regionAEnd

    orientation=$(echo $coords | awk '{print $3}')
    echo $orientation

    if [[ $orientation == '-' ]]; then
        echo ${file%.fa} >> /data/wrayva/output/extract_regions/regionA2/chm13chr4/reverseOrientedMatches.txt
    fi

    #write region A coords to bed file
    #echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName}
    echo -e $chrName"\t"$regionAStart"\t"$regionAEnd"\t"regionA_${hapName} > /data/wrayva/output/extract_regions/regionA2/chm13chr4/bed/regionA_${file%.fa}.bed

    #extract with samtools
    #cutStart="$(($regionAStart + 1))"
    #zeroIndexedStart="$(($regionAStart))"
    #cutEnd="$(($regionAEnd))"
    #zeroIndexedEnd=$cutEnd
    #echo ">${file%.fa}" > /data/wrayva/output/extract_regions/regionA2/chm13chr4/sequences/regA_${file}
    #cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/regionA2/chm13chr4/sequences/regA_${file}
    samtools faidx $file ${file%.fa}:${regionAStart}-${regionAEnd} > /data/wrayva/output/extract_regions/regionA2/chm13chr4/sequences/regA_${file}

    #output region A length to file
    length="$(($regionAEnd - $regionAStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/regionA2/chm13chr4/lengths.txt
done

#testing
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4/full_sequences
for file in `ls *.fa | sort -R | head -n1`; do
    echo ${file%.fa}
    genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
    chrName=$(echo ${file%.fa} | cut -d '_' -f2)
    hapName=$(echo ${file%.fa} | cut -d '_' -f3)
    echo "hapName: $hapName"

    cat /data/wrayva/output/extract_regions/regionA2/chm13chr4/filtered_minimap_chr4_regA.paf | grep ${genomeName}_${chrName}_${hapName} | sort -nk10 | cut -c-100

    coords=$(cat /data/wrayva/output/extract_regions/regionA2/chm13chr4/filtered_minimap_chr4_regA.paf | grep ${genomeName}_${chrName}_${hapName} | sort -nrk10 | head -n1 | awk '{print $3,$4,$5}')

    regionAStart=$(echo $coords | awk '{print $1}')
    echo $regionAStart

    regionAEnd=$(echo $coords | awk '{print $2}')
    echo $regionAEnd

    orientation=$(echo $coords | awk '{print $3}')
    echo $orientation
done

#for now, consider reverse oriented haps to be bad
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4/sequences
for file in `ls regA_*fa`; do
    file2=${file#regA_}
    file2=${file2%.fa}
    if [[ $(cat ../reverseOrientedMatches.txt | grep ${file2} | wc -l) -ge 1 ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/regionA2/chm13chr4/bad_regionA.txt
    else
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1  >> /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA.txt
    fi
done

#put all good ones in one fasta
for name in `cat /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA.txt`; do
    cat /data/wrayva/output/extract_regions/regionA2/chm13chr4/sequences/regA_${name}.fa >> /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA.fna
done

toReplace=">chr3_pat_hsa4:264157-320643"
echo $toReplace

replaceWith=">mPanPan1_chr3"
echo $replaceWith

sed -i "559085s/$toReplace/$replaceWith/g" /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_with_bonobo_mafft_aln.aln

#/data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA.fna contains the good chr4 region A


#combine the chr4 region A with the acros reg A
cat /data/wrayva/output/extract_regions/regionA2/good_regionA.fna /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA.fna > /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA_combined.fna


sed -i

#!/bin/bash
set -e
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK --retree 2 --maxiterate 0 /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_combined_with_bonobo.fna > /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_with_bonobo_mafft_aln.aln

module load iqtree

iqtree2 -T $SLURM_CPUS_PER_TASK -s /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_with_bonobo_mafft_aln.aln

sbatch --time=8:00:00 --mem=8g --cpus-per-task=5 /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/run_mafft_and_iqtree.sh
#31256623

sbatch --time=5:00:00 --mem=8g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_with_bonobo_mafft_aln.aln /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/ "-o mPanPan1_chr3 -redo"

#remove distal_NA18983_chr4_haplotype2-0000076_0-28518; it only contains gaps

seqkit grep -rvip "^distal_NA18983_chr4_haplotype2-0000076:0-28518" good_regionA_combined.fna > good_regionA_combined_removedGapSeq.fna
seqkit grep -rvip "^distal_NA18983_chr4_haplotype2-0000076:0-28518" good_regionA_combined_mafft_aln.aln > good_regionA_combined_mafft_aln_removedGapSeq.aln

#rename distal_genomename_chr_... to chr_genomename_... to make it clearer in ncbi msa viewer
grep -n -e "distal_NA[0-9]*_chr4" good_regionA_combined_mafft_aln_removedGapSeq.aln

cp good_regionA_combined_mafft_aln.aln good_regionA_combined_mafft_aln_copy.aln

for file in `grep ">" good_regionA_combined_mafft_aln_copy.aln`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f1 | cut -d '>' -f2) #orig file had -f3
    hapName=$(echo ${file%.fna} | cut -d '_' -f3 | cut -d ':' -f1) #orig file had -f4
    echo "hapName: $hapName"

    toReplace=distal_${genomeName}_${chrName}
    echo $toReplace

    replaceWith=${chrName}_${genomeName}
    echo $replaceWith

    lineNumber=$(grep -n -e "distal_${genomeName}_${chrName}_${hapName}" good_regionA_combined_mafft_aln_copy.aln | cut -d ':' -f1)
    echo $lineNumber

    replaceWith2=">${chrName}_${genomeName}_${hapName}"

    #sed -i "${lineNumber}s/$toReplace/$replaceWith/g" good_regionA_combined_mafft_aln_copy.aln
    sed -i "${lineNumber}s/*/$replaceWith2/g" good_regionA_combined_mafft_aln_copy.aln
done

grep ">" good_regionA_combined_mafft_aln_removedGapSeq.aln

cp /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA_combined_mafft_aln_chopped.aln.treefile /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA_combined_mafft_aln_chopped_copy.aln.treefile

cp good_regionA_with_bonobo_mafft_aln.aln.treefile good_regionA_with_bonobo_mafft_aln_copy.aln.treefile

#rename to include superpop and pop
cd /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=/data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_with_bonobo_mafft_aln_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/wrayva/output/extract_regions/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo

#create tree plots
module load R
R

library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/good_regionA_with_bonobo_mafft_aln_copy_rescale.aln.treefile")
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

ggsave("/data/wrayva/output/plots/regionATree_withChr4_bonobo_rescale5_color_by_chr.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

pops <- list(EUR=c(tree$tip.label[grep("EUR", tree$tip.label)]),
             EAS=c(tree$tip.label[grep("EAS", tree$tip.label)]),
             AMR=c(tree$tip.label[grep("AMR", tree$tip.label)]),
             AFR=c(tree$tip.label[grep("AFR", tree$tip.label)]),
             SAS=c(tree$tip.label[grep("SAS", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, pops,
                       group_name = "population")
ggtree(grouped_tree2, aes(color=population), size=3 layout='circular') +
  theme(legend.position=c(.95,.95), legend.key.size = unit(.2, 'cm')) +
  geom_tiplab(align=TRUE, linesize=.5, size=4, fontface='bold')
ggsave("/data/wrayva/output/plots/regionATree_withChr4_bonobo_rescale2_color_by_super_pop.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

#
#
#chop beginning of MSA off at 32690 to avoid building tree on MSA with lots of gaps
#
module load samtools
cd /data/wrayva/output/extract_regions/regionA2/chm13chr4
for file in `grep ">" good_regionA_combined_mafft_aln_copy.aln`; do
    echo ${file#>}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f1 | cut -d '>' -f2) #orig file had -f3
    hapName=$(echo ${file%.fna} | cut -d '_' -f3) #orig file had -f4
    #echo "hapName: $hapName"
    #echo $genomeName
    #echo $chrName

    samtools faidx good_regionA_combined_mafft_aln_copy.aln ${file#>}:32690-77584 >> good_regionA_combined_mafft_aln_chopped.aln
done

45750

#run iqtree
sbatch --time=8:00:00 --mem=16g --cpus-per-task=5 /data/wrayva/scripts/run_iqtree.sh /data/wrayva/output/extract_regions/regionA2/chm13chr4/good_regionA_combined_mafft_aln_chopped.aln /data/wrayva/output/extract_regions/regionA2/chm13chr4/
#31330144


#find bonobo haplotype
find -name \'\*bonobo*\'
find -name \'\*pan*\'
find . -name "*pan*"

here: (zipped)
/data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0.pri.fa.gz
or (not zipped)
/data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa

#align chr4regA to bonobo
#!/bin/bash
set -e
module load minimap2
minimap2 -cx asm20 /data/wrayva/output/extract_regions/regionA2/chm13chr4/chm13chr4regionA.fa /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa > minimap_chr4_regA_bonobo.paf

sbatch --time=5:00:00 --mem=32g /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/run_minimap_bonobo.sh
#31325523

#filter paf
awk '($17~"tp:A:P")  {print $0}' /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/minimap_chr4_regA_bonobo.paf > /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/filtered_minimap_chr4_regA_bonobo.paf

#extract bonobo reg A from minimap alignment coords
coords=$(cat /data/wrayva/output/extract_regions/regionA2/chm13chr4/bonobo/filtered_minimap_chr4_regA_bonobo.paf | sort -nrk10 | head -n1 | awk '{print $1,$3,$4,$5}')

sequenceName=$(echo $coords | awk '{print $1}')
echo $sequenceName

regionAStart=$(echo $coords | awk '{print $2}')
echo $regionAStart

regionAEnd=$(echo $coords | awk '{print $3}')
echo $regionAEnd

orientation=$(echo $coords | awk '{print $4}')
echo $orientation

samtools faidx /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa ${sequenceName}:${regionAStart}-${regionAEnd} >> good_regionA_combined_with_bonobo.fna

#extract region A chm13 chr22
samtools faidx /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa chr22:4578-67032 > /data/wrayva/output/extract_regions/regionA2/regionA.fa

#look at bad region A, see how many are actually missing region A vs just weird coordinates

module load samtools

origRegion=regionA
region=regionA2
extractDir=/data/wrayva/output/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
sequencesDir=/data/wrayva/output/sequences

cd $regionDir
for file in `cat /data/wrayva/output/extract_regions/regionA2/census/noRegionA.txt`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"

    cat ../../sequences/${file}.fna >> bad_regionA2.fna
done

#align all bad region A sequences to masked chr22 with minimap
sbatch --time=2:00:00 --mem=32g --cpus-per-task=20 ${scriptDir}/run_minimap.sh /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/distal_masked_refs/chr22_distal/chr22_distal.mask.fa bad_regionA2.fna $regionDir/paf/bad_regionA_minimap_aln2.paf
#31792754
#31882815

#filter minimap output
awk '($17~"tp:A:P")  {print $0}' $regionDir/paf/bad_regionA_minimap_aln2.paf > $regionDir/paf/filtered_bad_regionA_minimap_aln2.paf

#replace chr22_distal_mask with chr22 for viewing in IGV
sed -i "s/chr22_distal_mask/chr22/g" $regionDir/paf/filtered_bad_regionA_minimap_aln2.paf

#view the above paf in IGV

#moddotplot some of the bad ones
cd /data/wrayva/gitRepos/ModDotPlot
source venv/bin/activate
moddotplot interactive -f /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa /data/wrayva/output/sequences/distal_HG04199_chr13_haplotype1-0000025.fna -o /data/wrayva/output/moddotplot/interactive --compare-only --port $PORT1

moddotplot interactive -f /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa /data/wrayva/output/sequences/distal_NA18522_chr22_haplotype1-0000006.fna -o /data/wrayva/output/moddotplot/interactive --compare-only --port $PORT1

moddotplot interactive -f /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa /data/wrayva/output/sequences/distal_HG02559_chr21_haplotype2-0000222.fna -o /data/wrayva/output/moddotplot/interactive --compare-only --port $PORT1

#loop through sequences, check if properly assembled (contains telomere on one side, rdna on the other)
#check if contains region A
module load seqtk
cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"

    #check if contains region A
    foundRegionA=0
    teloPosition=$(seqtk telo $file | head -n1 | awk '{print $3}')
    if [[ -z "$teloPosition" ]]; then
        echo "no telomere"
        echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionA2/census/notProperlyAssembled.txt
    else
        echo "checking for region A"

        ${scriptDir}/run_minimap.sh $file $regionDir/chr22_regionA.fna $regionDir/paf/temp.paf

        awk '($17~"tp:A:P")  {print $0}' $regionDir/paf/temp.paf > $regionDir/paf/filtered_temp.paf

        #cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf

        #foundRegionA=0
        cat $regionDir/paf/filtered_temp.paf | while read -r line; do
            start=$(echo $line | awk '{print $3}')
            end=$(echo $line | awk '{print $4}')
            length=$((end - start))
            echo $length
            if [[ $length -ge 40000 ]]; then
                echo "found region A, length $length"
                #foundRegionA=1
                echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionA2/census/containsRegionA.txt
                break
            fi
        done
    fi
done

cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"

    if [[ -z $(grep ${file%.fna} /data/wrayva/output/extract_regions/regionA2/census/containsRegionA.txt) && -z  $(grep ${file%.fna} /data/wrayva/output/extract_regions/regionA2/census/notProperlyAssembled.txt) ]]; then
        echo "does not contain region A"
        echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionA2/census/noRegionA3.txt
    fi
done

for file in `cat /data/wrayva/output/extract_regions/regionA2/census/containsRegionA.txt`; do
    #echo $file
    #grep $file $regionDir/paf/filtered_bad_regionA_minimap_aln2.paf >> $regionDir/paf/filtered_bad_regionA_minimap_aln3.paf
    grep $file /data/wrayva/output/extract_regions/regionA2/census/notProperlyAssembled.txt
done

sequencesDir=/data/wrayva/output/sequences
cd $sequencesDir
ls distal_*fna > /data/wrayva/output/extract_regions/regionA2/census/total_assembled.txt

genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
outputDir=/data/wrayva/output/extract_regions/regionA2/census
cd $genomeDir
#loop through genomes
echo "SampleID, Number of complete distal regions identified, Number of distal regions containing A" > $outputDir/regionACensusByIndividual.csv
count=0
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    #echo $genomeName
    totalComplete=$(grep -c $genomeName $outputDir/total_assembled.txt)
    containingA=$(grep -c $genomeName $outputDir/containsRegionA.txt)
    echo "${genomeName},${totalComplete},${containingA}" >> $outputDir/regionACensusByIndividual.csv
    ((count=$count + $containingA))
    echo $count
done
echo $count

# check identity percentage with chr4 region A
${scriptDir}/run_minimap.sh /data/wrayva/output/extract_regions/regionA2/chr22_regionA.fna /data/wrayva/output/extract_regions/regionA2/chm13chr4/chm13chr4regionA.fa /data/wrayva/output/extract_regions/regionA2/chm13chr4/minimapWithAcroRegAToGetIdentityPercentage.paf
