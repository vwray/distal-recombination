This document contains details on the analysis I have run as part of this project.

## Trimming telomere
The below command runs over fasta files in a directory, finds the position where the telomere ends on each sequence, and trims the sequence from this position.
```sh
module load seqtk
for file in `ls rdna_trim_distal_*fna`; do
    #Find the position of telomere
    seqtk telo $file > /data/wrayva/output/seqtk_telo/output_${file%.fna}.txt
    teloPosition=$(cat /data/wrayva/output/seqtk_telo/output_${file%.fna}.txt | head -n1 | awk -F$'\t' '{print $3}')
    echo $teloPosition
    #trim the telomere from the sequence
    seqtk trimfq -b $teloPosition $file > /data/wrayva/output/trim_telo/trim_telo_${file%.fna}.fna
done
```

## Aligning Distal Regions
### parsnp
The below command is used to align each distal bit fasta sequence to the distal bit of chromosome 15 in CHM13 using `parsnp`.
```sh
parsnp -r /data/Phillippy2/projects/acro_comparisons/refs/CHM13/distal_bits/chr15.distal.fa -d /data/wrayva/output/sequences/ -o /data/wrayva/output/parsnp-out
```

### wfmash

The below command is used to align each distal bit fasta sequence to the distal bit of chromosome 15 in CHM13 using `wfmash`.
```sh
module load wfmash
wfmash --nosplit -s 100k /data/wrayva/output/chm13ref/trim_telo_trim_rdna_chr15.distal.fa /data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa > /data/wrayva/output/wfmash_on_trimmed_to_chm13/aln.paf
```

### Filter extracted DJ
The below bash code filters out extracted DJ sequences that do not match the first and last 10 characters of the expected DJ.
```sh
#good and bad DJ
for file in `ls /data/wrayva/output/extract_regions/DJ/sequences/dj_distal_*fna`; do
    if [[ $(cat $file | head -n2 | tail -n1 | cut -c-10) == 'AATGTAAGGC' && $(cat $file | head -n2 | tail -n1 | rev | cut -c-10) == 'CGGAAGTTCA' ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/DJ/good_DJ.txt
    else
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> /data/wrayva/output/extract_regions/DJ/bad_DJ.txt
    fi
done
```


# Fasta Files
## Rename Sequence Identifiers in Fasta File:
I used the following bash code to replace sequence identifiers like sequence_100 with distal_genomeName_chrName_hapName.
```bash
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
```

## Copy All Sequences to One Fasta File:
I used the following bash code to copy all sequences in `/data/wrayva/output/extract_regions/DJ` to one fasta file.
```bash
cd /data/wrayva/output/extract_regions/DJ
for file in `ls dj_distal_*fna`; do
    cat $file >> /data/wrayva/output/extract_regions/DJ/all_dj.fna
done
```

## Convert FASTA file to upper case:
The below bash code converts the specified FASTA file to uppercase (including the sequence identifiers):
```bash
awk '{ print toupper($0) }' /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/distal_masked_refs/chr22_distal/chr22_distal.mask.fa > /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.mask.upper.fa
```

## Index a FASTA file:
The below bash code uses samtools to index a fasta file.
```bash
module load samtools
samtools faidx /data/wrayva/output/sequences/all.fna distal_HG02273_chr22_haplotype1-0000003 > distal_HG02273_chr22_haplotype1-0000003.fa
```

## Get first 100 characters in FASTA file:
The below bash code gets the second line of a fasta file (where the sequence is) and extracts the first 100 characters.
```bash
cat /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa | head -n2 | tail -n1 | cut -c-100
```

## Get last 100 characters in FASTA file:
The below bash code gets the second line of a fasta file (where the sequence is) and extracts the last 100 characters (in reverse order).
```bash
cat /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa | head -n2 | tail -n1 | rev | cut -c-100
```

## Extract and reverse complement a region of a haplotype
The below bash code extracts the specified sequence at the specified coordinates from the specified file, doing the reverse complement and not marking the sequence with "\rc".
```bash
samtools faidx /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa ${sequenceName}:${bonoboStart}-${bonoboEnd} --reverse-complement --mark-strand no >> ${regionDir}/bonobo/bonobo.fna
```

# PAF Files
## Create paf with a selection of haplotypes:
The below bash code extracts certain haplotypes listed in text file to a new paf file for viewing in IGV.
```bash
cat /data/wrayva/output/extract_regions/DJ/too_small_dj.txt | while read -r line; do
    echo $line
    genomeName=$(echo $line | awk '{print $1}')
    chrName=$(echo $line | awk '{print $3}')
    hapName=$(echo $line | awk '{print $2}')
    echo "hapName: $hapName"
    cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | grep distal_${genomeName}_${chrName}_${hapName} >> /data/wrayva/output/extract_regions/DJ/validation/bad_djs.paf
done
```

## Replace a word/phrase in a file:
The below bash code replaces every occurrence of "chr22_distal_mask" with "chr22" in the specified file.

```bash
sed -i 's/chr22_distal_mask/chr22/g' /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf
```

## Filter a PAF file:
The below bash code filters a PAF file to filter out secondary matches; only include primary and supplementary. If the PAF file does not include cigar strings, column 17 should be changed to column 13.
```bash
awk '($17~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_chm13chr22masked/minimap_output_with_cigar.paf > /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf
```

## Query a PAF file to get overlapping regions:
The below bash code searches a PAF file for regions overlapping with the desired start and end, matching the desired haplotype name, sorted by haplotype, first 100 rows.
```bash
cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_aln.paf | awk '$8 <= 67033 && $9 >= 4577' | grep distal_HG00099_chr13_haplotype1-0000019 | sort -k1 | head -n100
```

## Look at PAF file in terminal cutting off cigar strings:
The below bash code displays PAF file contents in the terminal, with the junky columns 19 through 24 removed, and sorted numerically by column 8 (reference start coordinate).
```bash
cat chm13chr22vschr21.paf | sort -nk8 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}'
```

#SAM/BAM Files
## Filter BAM file based on region in a bed file:
The below bash code filters a BAM file, returning only alignments that overlap with the region in the bed file.
```bash
module load samtools
samtools view -L regionA.bed /data/wrayva/output/minimap_chm13chr22masked/filtered.bam
```

## Sort and convert a SAM file to a compressed BAM file
```bash
samtools sort  -T /data/wrayva/output/minimap_chm13chr22masked/sam_to_bam.tmp -O bam /data/wrayva/output/minimap_chm13chr22masked/minimap_output_as_sam.sam > /data/wrayva/output/minimap_chm13chr22masked/minimap_output_as_bam.bam
```

## Index a BAM file
```bash
samtools index /data/wrayva/output/minimap_chm13chr22masked/minimap_output_as_bam.bam
```

#Plots
##ModDotPlot
The below code uses a github installation of [ModDotPlot](https://github.com/marbl/ModDotPlot) with Conda to generate heatmaps for visualizing alignments of repetitive sequences.
```bash
source myconda
conda activate dfv
cd /data/wrayva/gitRepos/ModDotPlot
moddotplot static -f /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa /data/wrayva/output/sequences/distal_HG03831_chr14_haplotype1-0000003.fna -o /data/wrayva/output/moddotplot/HG03831_chr14_haplotype1-0000003 --compare
```

##Python

### Create histogram of distal bit lengths:
This bash code first creates a file `distalbits_lineLengths.txt` containing the line lengths of each line in the fasta file. Then we remove the line lengths corresponding to sequence identifier headings to create `distalbits_modifiedlineLengths.txt`.
```bash
awk '{ print length }' distalbits.fasta >> distalbits_lineLengths.txt

i=1
for line in `cat distalbits_lineLengths.txt`; do
    if [[ $(($i % 2)) -eq 0 ]]; then
        echo $line >> distalbits_modifiedlineLengths.txt
    fi
    ((i++))
done
```
This python code creates a histogram from the line lengths, which are the lengths of the distal regions:
```python
import numpy as np
import matplotlib.pyplot as plt
file_obj = open("/path/distalbits_modifiedlineLengths.txt", "r")

# create a list containing each length from the file, which are all on separate lines in the file
file_data = file_obj.read()
lines = file_data.splitlines()
x=[]
for i in range(0, len(lines)):
    x.append(int(lines[i]))

fig, ax = plt.subplots(1, 1)
ax.hist(x)

# Add title and axis labels
ax.set_title("Distribution of Lengths of Distal Bits")=
ax.set_xlabel('Length of distal bit')
ax.set_ylabel('Frequency')

plt.show()
plt.savefig('/path/distalbits_lengthsHistogram.png')
```

### Create %GC Content Plot:
This Python script creates a plot of the %GC content of a random selection of 40 of the sequences.
```python
from Bio.SeqUtils import gc_fraction
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random

def compute_gc_content(record, window_size, step_size, ax, offset):
    """Plot windowed GC content on a designated Matplotlib ax."""
    #Computes GC content based on window and step size
    #Add the offset in order to stack sequence plots
    y = [100.0*gc_fraction(record[i : i + window_size])+offset for i in range(0,len(record) - window_size, step_size)]
    x = [i for i in range(0,len(record) - window_size, step_size)]
    y2 = [offset for i in range(0,len(record) - window_size, step_size)]
    ax.fill_between(x, y, y2, alpha=0.5)

#Opens the fasta file
fastaObject = pysam.FastaFile('/data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa')

#Set up the title and axis labels
fig, ax = plt.subplots(1, 1)
ax.set_ylim(bottom=0,top=4000)
ax.set_ylabel("GC(%)")
ax.set_xlabel("Distal Region Position")
plt.yticks([])
plt.title("%GC Content for 40 Random Distal Regions\n(Window Size 1000, Step Size 100)")

#Chooses 40 random sequences in the fasta file
random_indices = set()
while(len(random_indices)< 40):
    random_indices.add(random.randint(0,len(fastaObject.references)-1))
random_indices_list = list(random_indices)

max_length = 0
for i in range(len(fastaObject.references)):
    j=random_indices_list[i]
    region = fastaObject.fetch(region=fastaObject.references[j])
    #Add the sequence name label to the plot
    ax.annotate(text=fastaObject.references[j][7:24]+fastaObject.references[j][30], xy=(-1400000,i*100+10), fontsize=0.5)
    compute_gc_content(region, 1000, 100, ax, i*100)
    if len(region) > max_length:
        max_length = len(region)

ax.set_xlim([-1500000,max_length+100000])
fig.savefig('/data/wrayva/output/plots/gc_plot_1.pdf')
```

## R
### Create a circular phylogenetic tree
#### Colored by chromosome
```r
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft.fna.treefile")
pdf(file = "/data/wrayva/output/plots/regionATree_circular_color_branch_lengths2.pdf")
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
```

#### Colored by population
```r
#color by population
library(tidyverse)
library(ggtree)
tree <- read.tree("/data/wrayva/output/extract_regions/regionA/good_regionA_alignment_mafft_COPY.fna.treefile")
pdf(file = "/data/wrayva/output/plots/regionATree_circular_color_by_pop3.pdf")
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
```
