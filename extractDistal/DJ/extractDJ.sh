module load minimap2

cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"
    #get telomere coords from seqtk telo
    teloPosition=$(cat /data/wrayva/output/seqtk_telo/output_rdna_trim_${file%.fna}.txt | head -n1 | awk -F$'\t' '{print $3}')
    #handle empty minimap output
    teloFile="/data/wrayva/output/seqtk_telo/output_rdna_trim_${file%.fna}.txt"; [ -e $teloFile ] && [ -s $teloFile ] || teloPosition=0
    #write telomere coords to bed file
    echo -e $chrName"\t"0"\t"$teloPosition"\t"telo_${hapName} >> /data/wrayva/output/bed/${genomeName}.bed

    #run minimap on all sequences with my seq as ref and the DJ as query.
    minimap2 -x asm20 $file /data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_5S_rDNA_PHR/DJ.fa > /data/wrayva/output/minimap_dj/${file%.fna}.paf
    djStart=$(awk '($13~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_dj/${file%.fna}.paf | awk '$12 >= 60'| sort -nrk8 | head -n1 | awk '{print $8}')
    #grab entry 9, tab separated, in the paf file
    djEnd=$(awk '($13~"tp:A:P")  {print $0}' /data/wrayva/output/minimap_dj/${file%.fna}.paf | awk '$12 >= 60'| sort -nrk8 | head -n1 | awk '{print $9}')
    #write DJ coords to bed file
    echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName} >> /data/wrayva/output/bed/${genomeName}.bed
    # get lengths
    length="$((djEnd - teloPosition))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/bed/lengths.txt

    djBedLine=$(cat /data/wrayva/output/bed/${genomeName}.bed | grep DJ_${hapName})
    djStart=$(echo $djBedLine | awk '{print $2}')
    djEnd=$(echo $djBedLine | awk '{print $3}')

    #extract DJ
    #write DJ coords to new bedfile
    echo $djBedLine >> /data/wrayva/output/extract_regions/DJ/${file%.fna}.bed

    #extract with seqtk
    #seqtk subseq $file ${file%.fna}.bed > /data/wrayva/output/extract_regions/DJ/dj_${file}

    #extract with cut
    cutStart="$(($djStart - 1))"
    zeroIndexedStart="$(($cutStart - 1))"
    cutEnd="$(($djEnd + 6))"
    zeroIndexedEnd=$cutEnd
    echo ">${file%.fna}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/DJ/dj_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/DJ/dj_${file}

    #output DJ length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/DJ/lengths.txt
done
