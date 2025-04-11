This document contains documentation on the analysis I have run as part of this project using the V4 assemblies, as well as directories of results and data.

## Directories
#### Assemblies (Genome Directory)
`/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample`

#### Assemblies Graphs 
`/vf/users/Phillippy/projects/hprc-assemblies/assemblies-v4/HG00126/verkko-hi-c/`

#### Fasta Files
`/data/Phillippy2/projects/acro_comparisons/hprc/distal/v4/distalBitSequences/complete` 
This directory contains the fasta files of the complete distal sequences from the v4 assemblies, with telomere, rdna, and gap coordinates. Complete sequences are the ones that have telomere and rdna (more details on complete sequences below). Each file contains the first 7Mb of each haplotype.

#### Spreadsheet
`/data/Phillippy2/projects/acro_comparisons/hprc/distal/v4/distalBitSequences`
This directory contains the spreadsheet with coordinates of telomere, rDNA, and gaps as detailed below: `v4_assemblies_spreadsheet2.csv`
The filtered version of the spreadsheet, only containing complete distal sequences, is here: `v4_assemblies_spreadsheet_filtered.csv`

## Determining Complete Distal Regions
Complete distal regions are determined as ones that pass the following checks:
1. Contains telomere, as determined by using `seqtk telo`. Initially was checking for an entry in `$genomeDir/$genomeName/verkko*/analysis/${genomeName}.telomere.bed` for the given haplotype which starts at 0, but these coordinates were not all accurate as some were referring to coordinates before the reorienting of the haplotype.
2. Contains rDNA, as determined by running minimap2 with the distal region sample as query and rotated reference rDNA as reference. In most scenarios, the haplotype will contain a gap of 'N' characters in place of the rDNA, with a small amount of rDNA before and after the gap.

## Filling out Spreadsheet
This is the script located in this Github repo used to fill out the spreadsheet: `scripts/fillOutSpreadsheet.sh`

All coordinates are with respect to 0 being the beginning of the p-arm. As stated above, complete sequences are the ones that have telomere and rdna. The spreadsheet only contains complete sequences; others have been excluded.

The spreadsheet has the following columns:

#### DistalSample
This is the name of the distal sample, which includes genome name, chromosome, and haplotype in a format like {genomeName}_chr{chromosome number}_haplotype{haplotype number}-{sequence identifier number from assembly file}.

#### TelomereStart
This is the start position of the telomere found in the distal sample. This should always be 0, even in the case where no telomere is found.

#### TelomereEnd
This is the end position of the telomere found in the distal sample. If no telomere is found, we set this to 0, and we plan to exclude such samples from further analysis.

#### rDnaBeforeGapStart
This is the beginning position of the first segment of rDNA found by running minimap2 with the distal region sample as query and rotated reference rDNA (KY962518-ROT) as reference. Ideally this will be the beginning of the rDNA which is found immediately preceding a gap of 'N' characters.

#### rDnaBeforeGapEnd
This is the end position of the first segment of rDNA found by running minimap2 with the distal region sample as query and rotated reference rDNA (KY962518-ROT) as reference. Ideally this will coincide with the beginning of a gap of 'N' characters.

#### gapStart
This is the beginning position of the first segment of gap ('NNNNNNN'). In some case, multiple gaps are found which are close together. In this case, this is the beginning position of the first gap found. Ideally this will coincide with the end position of the rDNA before the gap.

#### gapEnd
This is the end position of the last segment of gap ('NNNNNNN'). In some case, multiple gaps are found which are close together. In this case, this is the end position of the last gap found. Ideally this will coincide with the beginning position of the rDNA after the gap.

#### rDnaAfterGapStart
This is the beginning position of the last segment of rDNA found by running minimap2 with the distal region sample as query and rotated reference rDNA (KY962518-ROT) as reference. Ideally this will coincide with the end of a gap of 'N' characters.

#### rDnaAfterGapEnd
This is the end position of the last segment of rDNA found by running minimap2 with the distal region sample as query and rotated reference rDNA (KY962518-ROT) as reference. Ideally this will be the end of the rDNA which is found immediately after a gap of 'N' characters.

## To do
* Use spreadsheet to compute basic statistics, i.e. how long is DJ, how many rDNA arrays complete vs gap, how many rDNAs span gap, how many rDNAs rotated
* Check if distal sequences without gap are assembled across rDNA
* Rerun analysis which has been done for V3 and needs to be repeated for V4