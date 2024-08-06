

module load samtools
module load R

origRegion=regionA
region=regionA2
extractDir=/data/wrayva/output/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
sequencesDir=/data/wrayva/output/sequences

mkdir -p $regionDir $regionDir/bed $regionDir/sequences $regionDir/paf $regionDir/filtered_paf ${regionDir}/mafft $regionDir/plots/

extractStart=$(cat ${extractDir}/${origRegion}.bed | awk '{print $2}')
extractEnd=$(cat ${extractDir}/${origRegion}.bed | awk '{print $3}')
expectedLength=$((extractEnd-extractStart))

#extract region from chm13 chr22
samtools faidx /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa chr22:${extractStart}-${extractEnd} > $regionDir/chm13chr22${region}.fa

#extract sequences based on alignment
cd $sequencesDir
for file in `ls distal_*fna`; do
    echo ${file%.fna}
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    #echo "hapName: $hapName"

    #align all sequences to extracted region with minimap
    #use extracted chm13 chr22 region as query, current sequence as reference
    ${scriptDir}/run_minimap.sh $file $regionDir/chm13chr22${region}.fa $regionDir/paf/${file%.fna}_minimap_${region}_aln.paf

    #filter minimap output
    awk '($17~"tp:A:P")  {print $0}' $regionDir/paf/${file%.fna}_minimap_${region}_aln.paf > $regionDir/filtered_paf/${file%.fna}_filtered_minimap_${region}_aln.paf

    coords=$(cat $regionDir/${file%.fna}_filtered_minimap_${region}_aln.paf | grep ${genomeName}_${chrName}_${hapName} | sort -nk3 | head -n1 | awk '{print $5,$8}')

    regionStart=$(echo $coords | awk '{print $2}')
    echo $regionStart

    orientation1=$(echo $coords | awk '{print $1}')

    coords=$(cat $regionDir/${file%.fna}_filtered_minimap_${region}_aln.paf | grep ${genomeName}_${chrName}_${hapName} | sort -nrk4 | head -n1 | awk '{print $5,$9}')

    regionEnd=$(echo $coords | awk '{print $2}')
    echo $regionEnd

    orientation2=$(echo $coords | awk '{print $1}')

    if [[ $orientation1 == '-' || $orientation2 == '-' ]]; then
        echo ${file%.fa} >> ${regionDir}/reverseOrientedMatches.txt
    fi

    #write coords to bed file
    echo -e $chrName"\t"$regionStart"\t"$regionEnd"\t"${region}_${hapName} > ${regionDir}/bed/${region}_${file%.fna}.bed

    #extract with samtools
    samtools faidx $file ${file%.fna}:${regionStart}-${regionEnd} > ${regionDir}/sequences/${region}_${file}

    #output region length to file
    length="$(($regionEnd - $regionStart))"
    echo -e $genomeName"\t"$hapName"\t"$chrName"\t"$length >> ${regionDir}/lengths.txt
done

#take "good" haplotypes
cd ${regionDir}/sequences
for file in `ls ${region}_*fna`; do
    genomeName=$(echo ${file%.fna} | cut -d '_' -f3)
    chrName=$(echo ${file%.fna} | cut -d '_' -f4)
    hapName=$(echo ${file%.fna} | cut -d '_' -f5)

    file2=${file#${region}_}
    file2=${file2%.fna}
    length=$(cat ../lengths.txt | grep $genomeName | grep $hapName | grep $chrName | awk '{print $4}')

    if [[ $(cat ../reverseOrientedMatches.txt | grep ${file2} | wc -l) -ge 1 ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> $regionDir/bad_${region}2.txt
    #figure out how to not hardcode the length cutoffs
elif [[ $length -le $((expectedLength - 10000)) || $length -ge $((expectedLength + 10000)) ]]; then
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1 >> $regionDir/bad_${region}2.txt
    else
        cat $file | head -n1 | cut -c2- | cut -d ':' -f1  >> $regionDir/good_${region}2.txt
    fi
done

#290 good ones

#put all good ones in one fasta
cd ${regionDir}
for file in `cat good_${region}.txt`; do
    cat $sequencesDir/${file}.fna >> ${regionDir}/sequences/all_${region}.fna
done

#fasta file with good sequences: ${regionDir}/sequences/all_${region}.fna

#check for bonobo region B
sbatch --time=5:00:00 --mem=32g --cpus-per-task=5 ${scriptDir}/run_minimap.sh /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa $regionDir/chm13chr22${region}.fa $regionDir/paf/bonobo_minimap_${region}_aln.paf
#31658209

#filter paf
awk '($17~"tp:A:P")  {print $0}' $regionDir/paf/bonobo_minimap_${region}_aln.paf > $regionDir/paf/filtered_bonobo_minimap_${region}_aln.paf

#extract bonobo region B to include in tree, add to ${regionDir}/sequences/all_${region}.fna
cat $regionDir/paf/bonobo_minimap_${region}_aln.paf
#chr22:226744-459516	232773	5	232769	+	chr14_pat_hsa13	116126430	8749771	8976845	139590

bonobo_coords=$(cat $regionDir/paf/filtered_bonobo_minimap_${region}_aln.paf | sort -nrk10 | head -n1 | awk '{print $6,$8,$9,$5}')

sequenceName=$(echo $bonobo_coords | awk '{print $1}')
echo $sequenceName

bonoboStart=$(echo $bonobo_coords | awk '{print $2}')
echo $bonoboStart

bonoboEnd=$(echo $bonobo_coords | awk '{print $3}')
echo $bonoboEnd

orientation=$(echo $bonobo_coords | awk '{print $4}')
echo $orientation

samtools faidx /data/Phillippy2/projects/primate_T2T/polishing/assemblies/mPanPan1_v2.0/mPanPan1_v2.0_pri.fa ${sequenceName}:${bonoboStart}-${bonoboEnd} >> ${regionDir}/sequences/all_${region}.fna

#run mafft on ${regionDir}/sequences/all_${region}.fna
cd ${regionDir}/mafft
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK --retree 2 --maxiterate 0 ${regionDir}/sequences/all_${region}.fna > ${regionDir}/mafft/mafft_${region}.aln

toReplace=$(grep $sequenceName ${regionDir}/sequences/all_${region}.fna)
echo $toReplace

replaceWith=">mPanPan1_$sequenceName"
echo $replaceWith

sed -i "s/$toReplace/$replaceWith/g" ${regionDir}/mafft/mafft_${region}.aln

module load iqtree
/data/wrayva/scripts/run_iqtree.sh ${regionDir}/mafft/mafft_${region}.aln ${regionDir}/mafft/ "-o $sequenceName"


#sbatch --time=5:00:00 --mem=32g --cpus-per-task=20 $regionDir/run_mafft_and_iqtree.sh $region $regionDir $sequenceName
#31709914

#rename things to have superpop and pop in name
cp ${regionDir}/mafft/mafft_${region}.aln.treefile ${regionDir}/mafft/mafft_${region}_copy.aln.treefile
cd ${genomeDir}
for genomeName in `ls -d HG* NA* | grep -v ".tar"`; do
    file=${regionDir}/mafft/mafft_${region}_copy.aln.treefile
    toReplace=distal_$genomeName

    pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    echo $superpop

    replaceWith=${superpop}_${pop}_${genomeName}

    sed -i "s/$toReplace/$replaceWith/g" $file
done

#need to resize bonobo branch
Rscript $scriptDir/plotTrees.r > $regionDir/plots/${region}_big_tree.pdf


#sbatch --time=5:00:00 --mem=32g --cpus-per-task=20 $scriptDir/alignAndCreateTrees50Random.sh
#31727316
#sbatch --time=2:00:00 --mem=16g --cpus-per-task=20 /data/wrayva/scripts/run_iqtree.sh ${subsetDir}/mafft/mafft_${region}.aln ${subsetDir}/mafft/ "-o $(echo $replaceWith | cut -c2-)"
#31731239

#run Parsnp
#sbatch --time=10:00:00 --mem=32g /data/wrayva/scripts/run_parsnp.sh /data/wrayva/output/extract_regions/regionA2/chr22_regionA.fna /data/wrayva/output/extract_regions/regionA2/good_sequences /data/wrayva/output/extract_regions/regionA2/parsnp-output -c
