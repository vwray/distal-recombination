cd /data/wrayva/output/sequences
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo "hapName: $hapName"


    #cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk8 | head -n1

    regionAStart=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577 && $3 < 100000' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nk3 | head -n1 | awk '{print $3}')
    echo $regionAStart

    regionAEnd=$(cat /data/wrayva/output/minimap_chm13chr22masked/filtered_minimap_with_cigar.paf | awk '$8 <= 67033 && $9 >= 4577 && $3 < 100000' | grep distal_${genomeName}_${chrName}_${hapName} | sort -nrk4 | head -n1 | awk '{print $4}')
    echo $regionAEnd


    if [[ $orientation1 == '-' || $orientation2 == '-' ]]; then
        echo ${file%.fna} >> /data/wrayva/output/extract_regions/regionA/reverseOrientedMatches.txt
    fi

    #write region A coords to bed file
    #echo -e $chrName"\t"$djStart"\t"$djEnd"\t"DJ_${hapName}
    echo -e $chrName"\t"$regionAStart"\t"$regionAEnd"\t"regionA_${hapName} > /data/wrayva/output/extract_regions/regionA/regionA_${file%.fna}.bed

    #extract with cut
    cutStart="$(($regionAStart + 1))"
    zeroIndexedStart="$(($regionAStart))"
    cutEnd="$(($regionAEnd))"
    zeroIndexedEnd=$cutEnd
    echo ">${file%.fna}:${zeroIndexedStart}-${zeroIndexedEnd}" > /data/wrayva/output/extract_regions/regionA/sequences/regA_${file}
    cat $file | head -n2 | tail -n1 | cut -c${cutStart}-${cutEnd} >> /data/wrayva/output/extract_regions/regionA/sequences/regA_${file}

    #output region A length to file
    length="$(($zeroIndexedEnd - $zeroIndexedStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> /data/wrayva/output/extract_regions/regionA/lengths.txt
done
