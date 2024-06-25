#!/bin/bash
#first argument is the directory where genomes are located
genomeDir=$1
#second argument is the output directory
outputDir=$2
#third argument is the script directory
scriptDir=$3

cd $genomeDir
#loop through genomes
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
  echo $genomeName
  cd $genomeDir/$genomeName/verkko-hi-c
  #get haplotypes that contain rdna
  for line in `grep 'rdna' assembly.paths.tsv | awk -F"\t" '{print $1}'`; do
    #get haplotype name from scfmap file
    hapName=$(grep -w $line assembly.scfmap | cut -d ' ' -f2)
    echo $hapName
    i=1
    #look up each haplotype name in the mashmap file, filtering on 99% coverage and length of match of 200,000
    grep -w $hapName analysis/assembly-ref.norm.mashmap |sed s/id:f://g |awk '{if ($(NF-1) > 0.99 && $4-$3 > 200000) print $0}' | sort -nk3 | while read -r mashmapLine; do
      if [[ $i -eq 1 && `echo $mashmapLine | cut -d ' ' -f5` == '-' ]]; then
        #get the chromosome label and orientation from the q arm of the mashmap results
        chrLabel=$(echo $mashmapLine | cut -d ' ' -f6)
        orientation='-'
        #find hapName in assembly file
        python $scriptDir/find_distal_bits.py -f "$genomeDir/$genomeName/verkko-hi-c/assembly.fasta" -a "$hapName" -r 'True' -o "$outputDir/distal_${genomeName}_${chrLabel}_${hapName}.txt"
        # add to fasta file
        echo ">${genomeName}_${chrLabel}_${hapName}" >> $outputDir/distalbits.fasta
        cat "$outputDir/distal_${genomeName}_${chrLabel}_${hapName}.txt" >> $outputDir/distalbits.fasta
      fi
      if [[ `echo $mashmapLine | cut -d ' ' -f5` == '+' ]]; then
        echo "$genomeName $hapName contains a +"
        if [[ $i -eq 1 ]]; then
          echo "distal bit is at beginning here, need to handle appropriately"
          echo "skipping $genomeName $hapName"
          #get chromosome label from last row instead of first
          chrLabel=$(grep -w $hapName analysis/assembly-ref.norm.mashmap |sed s/id:f://g |awk '{if ($(NF-1) > 0.99 && $4-$3 > 200000) print $0}' | sort -k3 | tail -n1 | awk -F"\t" '{print $6}')
          #echo "chr label for $genomeName $hapName is: $chrLabel"
          #python /home/wrayva/scripts/find_distal_bits.py -f "/data/Phillippy2/projects/hprc-assemblies/assemblies-v3/$genomeName/verkko-hi-c/assembly.fasta" -a "$hapName" -r 'False' -o "/home/wrayva/output/distal_${genomeName}_${chrLabel}_${hapName}.txt"
          # add to fasta file
          #echo ">${genomeName}_${chrLabel}_${hapName}" >> /home/wrayva/output/distalbits.fasta
          #cat "$outputDir/distal_${genomeName}_${chrLabel}_${hapName}.txt" >> $outputDir/distalbits.fasta
        fi
      fi
      ((i++))
    done
  done
done
