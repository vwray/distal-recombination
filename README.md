This repo contains code to extract distal regions of human chromosomes, align them, and build phylogenetic trees.

# Extracting Distal Regions
Use the bash script [extract_distal_regions](extractDistal/extract_distal_regions.sh) to extract the distal regions of the human chromosomes. You must provide a directory where the genome assemblies are located, an output directory, and a script directory that contains the python script [find_distal_bits](extractDistal/find_distal_bits.py).

For example,
```bash
./extract_distal_regions.sh /path/genomeDirectory /path/outputDirectory /path/scriptDirectory
```

The script performs the following logic:
* For each of the genomes in this folder: /data/Phillippy2/projects/hprc-assemblies/assemblies-v3/ the haplotypes of interest are retrieved by finding all instances of 'rdna' in assembly.paths.tsv .
* The haplotype name for each of these haplotypes comes from the $genomeName/verkko-hi-c/assembly.scfmap file.
* Next, the haplotype name is found in the analysis/assembly-ref.norm.mashmap file, with some filtering provided by Serge to include only matches with above 99% coverage and of length of the match at least 200,000. The output looks something like this:
```
haplotype2-0000156	98725104	10000	17500000	-	chr15	99753195	82262810	99750900	39	17490000	30	0.998881	kc:f:1.15289
haplotype2-0000156	98725104	17730000	19370000	-	chr15	99753195	80382628	82032174	40	1649546	31	0.999109	kc:f:1.066
haplotype2-0000156	98725104	19540000	26590000	-	chr15	99753195	73159639	80212695	38	7053056	31	0.999222	kc:f:0.761286
haplotype2-0000156	98725104	26590000	27770000	-	chr15	99753195	71955069	73135412	40	1180343	31	0.99925	kc:f:0.986899
haplotype2-0000156	98725104	27830000	45960000	-	chr15	99753195	53765511	71895009	37	18130000	31	0.999158	kc:f:0.98252
haplotype2-0000156	98725104	45970000	69340000	-	chr15	99753195	30383977	53761641	40	23377664	31	0.999229	kc:f:0.834773
haplotype2-0000156	98725104	69460000	71080000	+	chr15	99753195	28333401	29954910	37	1621509	29	0.998794	kc:f:1.23567
haplotype2-0000156	98725104	71190000	72940000	-	chr15	99753195	26309957	28061881	39	1751924	30	0.999105	kc:f:1.10247
haplotype2-0000156	98725104	73080000	78250000	-	chr15	99753195	20967309	26135797	36	5170000	29	0.99886	kc:f:0.929478
haplotype2-0000156	98725104	78400000	79140000	-	chr15	99753195	20068804	20817149	37	748345	29	0.99869	kc:f:1.03551
haplotype2-0000156	98725104	80690000	80930000	+	chr15	99753195	19000655	19240179	40	240000	25	0.996741	kc:f:1.0969
haplotype2-0000156	98725104	81590000	81870000	+	chr15	99753195	19004315	19284450	27	280135	27	0.998007	kc:f:1.02132
haplotype2-0000156	98725104	82180000	82480000	-	chr15	99753195	17665714	17971062	37	305348	29	0.998717	kc:f:1.02248
haplotype2-0000156	98725104	84030000	84350000	-	chr15	99753195	15686940	16007089	38	320149	35	0.999698	kc:f:0.809221
haplotype2-0000156	98725104	84750000	85210000	-	chr15	99753195	14759318	15223704	38	464386	34	0.9996	kc:f:0.885866
haplotype2-0000156	98725104	85210000	85470000	-	chr15	99753195	14524419	14772100	39	260000	26	0.997552	kc:f:1.03881
haplotype2-0000156	98725104	86500000	87220000	-	chr15	99753195	13899649	14619387	40	720000	36	0.99972	kc:f:0.961757
haplotype2-0000156	98725104	90800000	91030000	-	chr15	99753195	8636366	8861643	34230000	27	0.998054	kc:f:0.352476
haplotype2-0000156	98725104	91420000	91660000	-	chr15	99753195	8011292	8249569	34240000	32	0.999335	kc:f:0.263762
haplotype2-0000156	98725104	92970000	93180000	-	chr15	99753195	6096637	6271052	40210000	29	0.998864	kc:f:0.566217
haplotype2-0000156	98725104	93940000	94170000	-	chr15	99753195	5143916	5374092	36230176	39	0.999866	kc:f:0.962141
haplotype2-0000156	98725104	94180000	94440000	-	chr15	99753195	4874142	5133964	40260000	40	0.999889	kc:f:0.840008
haplotype2-0000156	98725104	96170000	96510000	-	chr14	101161492	1309234	1649139	39340000	34	0.999598	kc:f:1.09141
```
*The rows are sorted in ascending order of the third column (the start positions), and we assume that the first row is the q arm if there is a '-' in the 5th column. We may switch this logic to look instead for the largest segment (difference between fourth and third columns). This is also the row we use to get the chromosome label (6th column).
*Once we have the chromosome label and orientation, we call the python script [find_distal_bits](extractDistal/find_distal_bits.py) to reverse complement a haplotype if necessary so that the distal end on the p arm is first, find the 1Mb gap that indicates the rDNA, and extract the region before the gap. Gaps are indicated in the fasta file with 'N' characters.
*Finally, we write each distal region to a fasta file with a sequence identifier containing the genome name, chromosome label, and haplotype name.

## Trimming rdna

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

# Aligning Distal Regions
## parsnp
The below command is used to align each distal bit fasta sequence to the distal bit of chromosome 15 in CHM13 using `parsnp`.
```sh
parsnp -r /data/Phillippy2/projects/acro_comparisons/refs/CHM13/distal_bits/chr15.distal.fa -d /data/wrayva/output/sequences/ -o /data/wrayva/output/parsnp-out
```

## wfmash

The below command is used to align each distal bit fasta sequence to the distal bit of chromosome 15 in CHM13 using `wfmash`.
```sh
module load wfmash
wfmash --nosplit -s 100k /data/wrayva/output/chm13ref/trim_telo_trim_rdna_chr15.distal.fa /data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa > /data/wrayva/output/wfmash_on_trimmed_to_chm13/aln.paf
```

# Building Phylogenetic Trees
