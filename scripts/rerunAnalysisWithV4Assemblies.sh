#!/bin/bash

module load minimap2
file='/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample/HG00097/verkko-hi-c/analysis/HG00097.refOriented.fasta.gz'
minimap2 -cx asm20 $file /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa > /data/wrayva/output/rerun_v4/HG00097_DJ.paf

#above is in /data/wrayva/gitRepos/distal-recombination/scripts/runMinimap.sh

sbatch --time=5:00:00 --mem=32g /data/wrayva/gitRepos/distal-recombination/scripts/runMinimap.sh

cd /data/wrayva/output/rerun_v4/extract
sbatch --time=12:00:00 --mem=16g --cpus-per-task=20 extract_distal.sh

# above runs out of memory

# try running the same that i ran before to extract distal bits
#!/bin/bash
#first argument is the directory where genomes are located
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
#second argument is the output directory
outputDir=/data/wrayva/output/rerun_v4/extract
#third argument is the script directory
scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts


samtools faidx $file chr22_haplotype1-0000001 > /data/wrayva/output/rerun_v4/extract/HG00097_chr22_hap1.fa


module load samtools

genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
cd $genomeDir
#loop through genomes
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
  echo $genomeName
  cd $genomeDir/$genomeName/verkko*/analysis
  #get haplotypes that contain rdna
  for haplotype in `cat ${genomeName}.ref.reorient.tsv | grep -E 'chr13|chr14|chr15|chr21|chr22' | grep -v "unassigned" | cut -d " " -f4`; do
    echo $haplotype
    samtools faidx ${genomeName}.refOriented.fasta.gz ${haplotype}:1-7000000 > /data/wrayva/output/rerun_v4/extract/${genomeName}_${haplotype}.fa
  done
done

#above has been run on all genomes in the genomeDir and extracted first 7 million base pairs of all acrocentric chromosomes to /data/wrayva/output/rerun_v4/extract

#check for telomere
#check if $genomeDir/$genomeName/verkko*/analysis/${genomeName}.telomere.bed contains an entry for the haplotype name (not including chromosome)
#loop over genome dir or extract dir?

genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
cd $genomeDir
#loop through genomes
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
  echo $genomeName
  cd $genomeDir/$genomeName/verkko*/analysis
  #get haplotypes that contain rdna
  for haplotype in `cat ${genomeName}.ref.reorient.tsv | grep -E 'chr13|chr14|chr15|chr21|chr22' | grep -v "unassigned" | cut -d " " -f4`; do
    echo $haplotype
    samtools faidx ${genomeName}.refOriented.fasta.gz ${haplotype}:1-7000000 > /data/wrayva/output/rerun_v4/extract/${genomeName}_${haplotype}.fa
  done
done

extractAnalysis=/data/wrayva/output/rerun_v4/extractAnalysis
extractDir=/data/wrayva/output/rerun_v4/extract
cd $extractDir
for file in `ls -d HG* NA*`; do
  #echo ${file%.fa}
  genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
  chrName=$(echo ${file%.fa} | cut -d '_' -f2)
  hapName=$(echo ${file%.fa} | cut -d '_' -f3)
  #echo $hapName
  hapWithTel=`cat $genomeDir/$genomeName/verkko*/analysis/${genomeName}.telomere.bed | grep $hapName | grep -w 0`
  #echo $hapWithTel
  trimmedHapWithTel=`echo $hapWithTel | head -n1 | cut -d " " -f1`
  if [ -n "${trimmedHapWithTel}" ]; then
    echo "${genomeName}_${chrName}_${trimmedHapWithTel} contains telomere"
    #echo ${genomeName}_${chrName}_${trimmedHapWithTel} >> $extractAnalysis/hapsWithTelomere.txt
  else
    echo "${genomeName}_${chrName}_${hapName} does NOT contain telomere"
    echo "${genomeName}_${chrName}_${hapName}" >> $extractAnalysis/hapsWithNoTelomere.txt
  fi
done
#sed -i '/^$/d' $extractAnalysis/hapsWithTelomere.txt

cat $extractAnalysis/hapsWithTelomere.txt

#get number of extractions
ls -d HG* NA* | wc -l
#2659

#get number of extractions with telomere
wc -l $extractAnalysis/hapsWithTelomere.txt
#1968

#get number of extractions with no telomere
wc -l $extractAnalysis/hapsWithNoTelomere.txt
#691

#check for rDNA by doing minimap alignment
#first check regular rDNA, then check rotated ref
regularRDNA=/vf/users/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/chm13.chr21_rdna_array.fa
rotatedRefRDNA=/vf/users/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/KY962518-ROT.fa

extractAnalysis=/data/wrayva/output/rerun_v4/extractAnalysis
extractDir=/data/wrayva/output/rerun_v4/extract
rdnaDir=${extractAnalysis}/rdna

module load minimap2
for file in `cat $extractAnalysis/hapsWithTelomere.txt | head -n3 | tail -n1`; do
  genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
  chrName=$(echo ${file%.fa} | cut -d '_' -f2)
  hapName=$(echo ${file%.fa} | cut -d '_' -f3)
  #run minimap to check for rDNA
  minimap2 -cx asm20 $regularRDNA "$extractDir/${file}.fa" > $rdnaDir/${genomeName}_${chrName}_${hapName}.paf
  #filter the paf file
  awk '($17~"tp:A:P")  {print $0}' $rdnaDir/${genomeName}_${chrName}_${hapName}.paf > $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}.paf
  #process the paf file

done

cat $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}.paf | sort -nk8 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}'

#switch query & ref
minimap2 -cx asm20 $regularRDNA "$extractDir/${file}.fa" > $rdnaDir/${genomeName}_${chrName}_${hapName}_2.paf
  #filter the paf file
awk '($17~"tp:A:P")  {print $0}' $rdnaDir/${genomeName}_${chrName}_${hapName}_2.paf > $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}_2.paf
cat $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}_2.paf | sort -nk8 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}'



#check for rdna gap
samtools faidx ${extractDir}/HG00097_chr13_haplotype1-0000018.fa "chr13_haplotype1-0000018:1-7000000":3614826-3615817


#check for telomere
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
file=HG00097_chr13_haplotype1-0000018.fa
genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
chrName=$(echo ${file%.fa} | cut -d '_' -f2)
hapName=$(echo ${file%.fa} | cut -d '_' -f3)

# extract distal (through PJ) of CHM13 chr 21
module load samtools
samtools faidx /data/Phillippy2/projects/acro_comparisons/refs/chm13/acros/chr21.fa chr21:1-7000000 > /data/wrayva/output/rerun_v4/chm13/chr32_distal.fa


#!/bin/bash
chm13dir=/data/wrayva/output/rerun_v4/chm13
file=HG00097_chr13_haplotype1-0000018.fa
genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
chrName=$(echo ${file%.fa} | cut -d '_' -f2)
hapName=$(echo ${file%.fa} | cut -d '_' -f3)
# align with minimap
module load minimap2
minimap2 -cx asm20 $chm13dir/chr21_distal.fa "$extractDir/${file}" > $chm13dir/${genomeName}_${chrName}_${hapName}_aligned_with_chm13chr21.paf
#filter the paf
awk '($17~"tp:A:P")  {print $0}' $chm13dir/${genomeName}_${chrName}_${hapName}_aligned_with_chm13chr21.paf > $chm13dir/filtered_${genomeName}_${chrName}_${hapName}_aligned_with_chm13chr21.paf

#run the above
sbatch --time=5:00:00 --mem=32g $chm13dir/runMinimap.sh
#47132140

# view paf in IGV
#replace chr21:1-7000000 with chr21
cp $chm13dir/filtered_${genomeName}_${chrName}_${hapName}_aligned_with_chm13chr21.paf $chm13dir/igv_filtered_${genomeName}_${chrName}_${hapName}_aligned_with_chm13chr21.paf
sed -i 's/chr21:1-7000000/chr21/g' $chm13dir/igv_filtered_${genomeName}_${chrName}_${hapName}_aligned_with_chm13chr21.paf

#GenBank ID for rDNA reference: KY962518.1

#use rotated ref
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
file=HG00097_chr13_haplotype1-0000018.fa
genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
chrName=$(echo ${file%.fa} | cut -d '_' -f2)
hapName=$(echo ${file%.fa} | cut -d '_' -f3)
minimap2 -cx asm20 $rotatedRefRDNA "$extractDir/${file}.fa" > $rdnaDir/rotated_${genomeName}_${chrName}_${hapName}.paf
#filter the paf file
awk '($17~"tp:A:P")  {print $0}' $rdnaDir/rotated_${genomeName}_${chrName}_${hapName}.paf > $rdnaDir/filtered_rotated_${genomeName}_${chrName}_${hapName}.paf
cat $rdnaDir/filtered_rotated_${genomeName}_${chrName}_${hapName}.paf | sort -nk3 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}'

#try the other way 
minimap2 -cx asm20 "$extractDir/${file}.fa" $rotatedRefRDNA > $rdnaDir/rotated_${genomeName}_${chrName}_${hapName}_2.paf
#filter the paf file
awk '($17~"tp:A:P")  {print $0}' $rdnaDir/rotated_${genomeName}_${chrName}_${hapName}_2.paf > $rdnaDir/filtered_rotated_${genomeName}_${chrName}_${hapName}_2.paf
cat $rdnaDir/filtered_rotated_${genomeName}_${chrName}_${hapName}_2.paf | sort -nk8 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}'

#creating spreadsheet
echo "DistalSample,TelomereStart,TelomereEnd,rDnaBeforeGapStart,rDnaBeforeGapEnd,GapStart,GapEnd,rDnaAfterGapStart,rDnaAfterGapEnd" >> v4_assemblies_spreadsheet.csv

#grab the first row - this should be the rdna before gap
line1=`cat $rdnaDir/filtered_rotated_${genomeName}_${chrName}_${hapName}.paf | sort -nk3 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}' | head -n1`
rDnaBeforeGapStart=`echo $line1 | cut -d ' ' -f3`
rDnaBeforeGapEnd=`echo $line1 | cut -d ' ' -f4`

#grab last row and call this rdna after gap
line1=`cat $rdnaDir/filtered_rotated_${genomeName}_${chrName}_${hapName}.paf | sort -nk3 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}' | tail -n1`
rDnaAfterGapStart=`echo $line1 | cut -d ' ' -f3`
rDnaAfterGapEnd=`echo $line1 | cut -d ' ' -f4`

#





#run below to fill out spreadsheet

scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts
sbatch --time=5:00:00 --mem=32g ${scriptDir}/fillOutSpreadsheet.sh

#copy/move to appropriate folder
targetDir=/data/Phillippy2/projects/acro_comparisons/hprc/distal/v4/distalBitSequences

chm13dir=/data/wrayva/output/rerun_v4/chm13
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
regularRDNA=/vf/users/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/chm13.chr21_rdna_array.fa
rotatedRefRDNA=/vf/users/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/KY962518-ROT.fa
extractAnalysis=/data/wrayva/output/rerun_v4/extractAnalysis
extractDir=/data/wrayva/output/rerun_v4/extract
rdnaDir=${extractAnalysis}/rdna
scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts

cd $extractAnalysis
for line in `cat v4_assemblies_spreadsheet.csv | tail -n +2 | head -n30`; do
  #columns in spreadsheet
  #${file},${telomereStart},${telomereEnd},${rDnaBeforeGapStart},${rDnaBeforeGapEnd},${gapStartEnd},${rDnaAfterGapStart},${rDnaAfterGapEnd},${mashmapRdnaStart},${mashmapRdnaEnd} 

  file=`echo $line | cut -d ',' -f1`
  #echo $file
  rDnaBeforeGapStart=`echo $line | cut -d ',' -f4`
  rDnaBeforeGapEnd=`echo $line | cut -d ',' -f5`
  gapStartEnd=`echo $line | cut -d ',' -f6`
  rDnaAfterGapStart=`echo $line | cut -d ',' -f7`
  rDnaAfterGapEnd=`echo $line | cut -d ',' -f8`
  mashmapRdnaStart=`echo $line | cut -d ',' -f9`
  mashmapRdnaEnd=`echo $line | cut -d ',' -f10`

  if [[ "$rDnaBeforeGapStart" == "0" ]]; then
    echo $file
  fi


done


file=HG00126_chr13_haplotype1-0000008

genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
chrName=$(echo ${file%.fa} | cut -d '_' -f2)
hapName=$(echo ${file%.fa} | cut -d '_' -f3)

cat $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}.paf | sort -nk8 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}'

file=HG00320_chr22_haplotype2-0000051
seqtk telo "$extractDir/${file}.fa"

file=HG03458_chr15_haplotype2-0000082
teloPosition=`seqtk telo "$extractDir/${file}.fa" | head -n1 | awk -F$'\t' '{print $3}'`
echo $teloPosition

if [ -z "${teloPosition}" ]; then
  echo "$teloPosition no telo"

  #echo ${genomeName}_${chrName}_${trimmedHapWithTel} >> $extractAnalysis/hapsWithTelomere.txt
else
  echo "$teloPosition has telo"
  #echo "${genomeName}_${chrName}_${hapName}" >> $extractAnalysis/hapsWithNoTelomere.txt
fi

#get the rows from spreadsheet that have telomere (do not have a 0 in 3rd coordinate)
cd $extractDir
awk '($3~"0")  {print $0}' v4_assemblies_spreadsheet2.csv

awk '{if ($(NF-1) > 0.99 && $4-$3 > 200000) print $0}'

for line in `cat v4_assemblies_spreadsheet2.csv`; do
  echo $line | awk -F',' '{if ($3 > 0 && $9 > 0) print $0}' | awk 'BEGIN{FS=OFS=","}{NF--; print}' | awk 'BEGIN{FS=OFS=","}{NF--; print}' >> v4_assemblies_spreadsheet_filtered.csv
done

for line in `cat v4_assemblies_spreadsheet_filtered.csv | tail -n +2`; do
  genomeName=$(echo ${line} | cut -d ',' -f1 | cut -d '_' -f1)
  chrName=$(echo ${line} | cut -d ',' -f1 | cut -d '_' -f2)
  hapName=$(echo ${line} | cut -d ',' -f1 | cut -d '_' -f3)
  file=$(echo ${line} | cut -d ',' -f1)
  cp "$extractDir/${file}.fa" /data/Phillippy2/projects/acro_comparisons/hprc/distal/v4/distalBitSequences/complete
done

#copy spreadsheets to shared location
cp "/data/wrayva/output/rerun_v4/extractAnalysis/v4_assemblies_spreadsheet.csv" /data/Phillippy2/projects/acro_comparisons/hprc/distal/v4/distalBitSequences
cp "/data/wrayva/output/rerun_v4/extractAnalysis/v4_assemblies_spreadsheet_filtered.csv" /data/Phillippy2/projects/acro_comparisons/hprc/distal/v4/distalBitSequences

