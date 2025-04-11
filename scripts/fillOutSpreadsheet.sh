#!/bin/bash
chm13dir=/data/wrayva/output/rerun_v4/chm13
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v4-best-per-sample
regularRDNA=/vf/users/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/chm13.chr21_rdna_array.fa
rotatedRefRDNA=/vf/users/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/KY962518-ROT.fa
extractAnalysis=/data/wrayva/output/rerun_v4/extractAnalysis
extractDir=/data/wrayva/output/rerun_v4/extract
rdnaDir=${extractAnalysis}/rdna
scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts
cd $extractAnalysis
module load minimap2
module load seqtk
echo DistalSample,TelomereStart,TelomereEnd,rDnaBeforeGapStart,rDnaBeforeGapEnd,GapStart,GapEnd,rDnaAfterGapStart,rDnaAfterGapEnd,mashmapRdnaStart,mashmapRdnaEnd > v4_assemblies_spreadsheet2.csv 
for file in `cat $extractAnalysis/allHaps.txt`; do

  genomeName=$(echo ${file%.fa} | cut -d '_' -f1)
  chrName=$(echo ${file%.fa} | cut -d '_' -f2)
  hapName=$(echo ${file%.fa} | cut -d '_' -f3)

  #get telomere start and end
  #hapWithTel=`cat $genomeDir/$genomeName/verkko*/analysis/${genomeName}.telomere.bed | grep $hapName | grep -w 0`
  #telomereStart=`echo $hapWithTel | awk '{print $2}'`
  #telomereEnd=`echo $hapWithTel | awk '{print $3}'`
  telomereStart=0
  #seqtk telo "$extractDir/${file}.fa" > /data/wrayva/output/seqtk_telo/output_${file%.fna}.txt
  telomereEnd=$(seqtk telo "$extractDir/${file}.fa" | head -n1 | awk -F$'\t' '{print $3}')
  #echo $teloPosition

  if [ -z "${telomereEnd}" ]; then
    telomereEnd=0
  fi

  #run minimap to check for rDNA
  #minimap2 ref.fa query.fq > approx-mapping.paf
  minimap2 -cx asm20 $rotatedRefRDNA "$extractDir/${file}.fa" > $rdnaDir/${genomeName}_${chrName}_${hapName}.paf
  #filter the paf file
  awk '($17~"tp:A:P")  {print $0}' $rdnaDir/${genomeName}_${chrName}_${hapName}.paf > $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}.paf
  #process the paf file
  #grab the first row - this should be the rdna before gap
  line1=`cat $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}.paf | sort -nk3 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}' | head -n1`
  rDnaBeforeGapStart=`echo $line1 | cut -d ' ' -f3`
  rDnaBeforeGapEnd=`echo $line1 | cut -d ' ' -f4`

  #grap last row and call this rdna after gap
  line1=`cat $rdnaDir/filtered_${genomeName}_${chrName}_${hapName}.paf | sort -nk3 | awk '{$19=$18=$20=$21=$22=$23=$24="";print $0}' | tail -n1`
  rDnaAfterGapStart=`echo $line1 | cut -d ' ' -f3`
  rDnaAfterGapEnd=`echo $line1 | cut -d ' ' -f4`

  #get gap
  #adapt python script extractDistal/find_distal_bits.py for this
  #replace sequence identifiers; remove the ":1-7000000"
  sequenceIdentifierLine=`cat "$extractDir/${file}.fa" | grep ">"`
  sequenceIdentifier="${sequenceIdentifierLine#>}"
  #newSequenceIdentifier="${oldSequenceIdentifier}:1-7000000"
  #sed -i "s/$oldSequenceIdentifier/$newSequenceIdentifier/g" "$extractDir/${file}.fa" 
  #python $scriptDir/check_rdna_gap.py -f "$extractDir/${file}.fa" -a "$sequenceIdentifier"

  #returns the beginning of the first gap and the end of the last gap
  gapStartEnd=`python $scriptDir/check_rdna_gap.py -f "$extractDir/${file}.fa" -a "$sequenceIdentifier"`

  #run Steven's mashmap script:
  bash /data/wrayva/gitRepos/distal-recombination/scripts/mashmap_45S_copy.sh "$extractDir/${file}.fa" $extractAnalysis/mashmap/$file
  mashmapRdnaStart=`cat $extractAnalysis/mashmap/${file}*.bed | awk '$5>=99' | sort -nk2 | head -n1 | awk 'BEGIN {OFS="\t"}; {print $2}'`
  mashmapRdnaEnd=`cat $extractAnalysis/mashmap/${file}*.bed | awk '$5>=99' | sort -nk3 | tail -n1 | awk 'BEGIN {OFS="\t"}; {print $3}'`

  echo ${file},${telomereStart},${telomereEnd},${rDnaBeforeGapStart},${rDnaBeforeGapEnd},${gapStartEnd},${rDnaAfterGapStart},${rDnaAfterGapEnd},${mashmapRdnaStart},${mashmapRdnaEnd} >> v4_assemblies_spreadsheet2.csv
done
