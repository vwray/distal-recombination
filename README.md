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
* The rows are sorted in ascending order of the third column (the start positions), and we assume that the first row is the q arm if there is a '-' in the 5th column. We may switch this logic to look instead for the largest segment (difference between fourth and third columns). This is also the row we use to get the chromosome label (6th column).
*Once we have the chromosome label and orientation, we call the python script [find_distal_bits](extractDistal/find_distal_bits.py) to reverse complement a haplotype if necessary so that the distal end on the p arm is first, find the 1Mb gap that indicates the rDNA, and extract the region before the gap. Gaps are indicated in the fasta file with 'N' characters.
* Finally, we write each distal region to a fasta file with a sequence identifier containing the genome name, chromosome label, and haplotype name.
* Note: I have verified that this chromosome assignment gives the same assignment as taking a consensus on the q-arm, as done in `analysis/reOrient*fasta`.

# Align distal bits to find intra-satellite segmental duplications
The following command is used to run minimap with the reference as CHM13 chromosome 22 with the satellites masked, and the query as a fasta file containing all the distal bit strings as extracted above, with potentially pieces of telomere and rDNA still intact on the ends. We are using chromosome 22 because it is the longest, so most likely to align with the most pieces. This overall alignment will allow us to find the positions of where each seg dup starts and ends in each haplotype in the query fasta file.
```sh
module load minimap2
minimap2 -cx asm20 /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/distal_masked_refs/chr22_distal/chr22_distal.mask.fa /data/wrayva/output/sequences/all.fna > /data/wrayva/output/minimap_chm13chr22masked/minimap_output_with_cigar.paf
```

The positions of the seg dups on CHM13 chromosome 22 can be found on the satellite annotation Cen/Sat v2.1 on the [Repeat annotation for CHM13](https://github.com/marbl/CHM13?tab=readme-ov-file#repeat-annotation). This satellite track can be pulled into the Integrative Genomics Viewer (IGV) by File -> Load from URL. The paf file generated from minimap can be pulled into IGV from File -> Load from File.

# Extract Segmental Duplications
We wish to exract four different regions of segmental duplications contained in the distal bits of the acrocentric chromosomes. Here are their coordinates in CHM13 chromosome 22, as given in the [satellite track](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed). Note that region D is the distal junction (DJ).

| Region | Start | End | Name |
| --- | --- | --- | --- |
| A | 4578 | 67032 | ct_22_1 |
| B | 349847 | 459516 | ct_22_5 |
| C | 3634820 | 3841139 | ct_22_6 through ct_22_10 |
| D | 4448073 | 4793794 | ct_22_15 through ct_22_18 |

Using the coordinates above to search the minimap PAF output, we can extract the regions A, B, C, and D as done in [extractDJ.sh](extractDistal/DJ/extractDJ.sh)

# Align Segmental Duplications
Once a particular region, say region A, is extracted from each haplotype, we can align these all to each other with [mafft](https://mafft.cbrc.jp/alignment/software/tips.html) to get a multiple sequence alignment (MSA). Use the script [run_mafft.sh](alignment/run_mafft.sh) for this. Since regions B and C are inverted repeats of each other, we will reverse complement every haplotype region matching region C (in the forward direction) and align these in one MSA together with region B.

# Building Phylogenetic Trees
Once the MSA is obtained, we can infer a phylogenetic tree using [IQ-TREE](http://www.iqtree.org/), which uses maximum likelihood. Gaps are treated as unknown characters, i.e. no information. Use the script [run_iqtree.sh](tree/run_iqtree.sh), passing in the alignment obtained from mafft.
