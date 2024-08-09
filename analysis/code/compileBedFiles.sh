#just do one big bed file with everything in it, to go with one fasta file containing all distal bits
#include population, super pop
#column 1 should match fasta.
#column 4 should include region _ super pop _ sub pop
#mention in slack message that the haplotype identifiers like haplotype2-0000116 may be different numbers than what they have.
region=DJ
extractDir=/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions
regionDir=${extractDir}/${region}
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
bedDir=/data/Phillippy2/t2t-share/distal-recombination/bed
sharedDir=/data/Phillippy2/t2t-share/distal-recombination

cd $regionDir

for file in `cat good_${region}.txt`; do
    genomeName=$(echo ${file} | cut -d '_' -f2)
    chrName=$(echo ${file} | cut -d '_' -f3)
    hapName=$(echo ${file} | cut -d '_' -f4)
    pop=$(cat /data/Phillippy2/projects/acro_comparisons/hprc/distal/populations.csv | grep $genomeName | cut -d ',' -f2)
    superpop=$(cat /data/nhansen/T2T_Globus_NFH_Archive/AnVIL_3202_samples/allele_frequencies/superpopulations.txt | grep $pop | awk -F"\t" '{print $7}')

    #look up bed file
    start=$(cat bed/*${file}.bed | awk '{print $2}')
    end=$(cat bed/*${file}.bed | awk '{print $3}')
    echo -e ${file}"\t"${start}"\t"${end}"\t"${region}_${pop}_${superpop} >> ${sharedDir}/hprc_distal_regions.bed
    #cat bed/*${file}.bed >> $bedDir/${genomeName}.bed
done

#check which haps I have all 4 regions for
sharedDir=/data/Phillippy2/t2t-share/distal-recombination
cd /data/Phillippy2/projects/acro_comparisons/hprc/distal/sequences
for file in `ls distal_*.fna`; do
    genomeName=$(echo ${file%.fna} | cut -d '_' -f2)
    chrName=$(echo ${file%.fna} | cut -d '_' -f3)
    hapName=$(echo ${file%.fna} | cut -d '_' -f4)
    echo $hapName

    regA=$(cat ${sharedDir}/hprc_distal_regions.bed | grep ${file%.fna} | grep regionA)
    echo $regA
    regB=$(cat ${sharedDir}/hprc_distal_regions.bed | grep ${file%.fna} | grep regionB)
    regC=$(cat ${sharedDir}/hprc_distal_regions.bed | grep ${file%.fna} | grep regionC)
    dj=$(cat ${sharedDir}/hprc_distal_regions.bed | grep ${file%.fna} | grep DJ)

    if [[ -n $regA && -n $regB && -n $regC && -n $dj ]]; then
        echo ${file%.fna} >> $sharedDir/hapsContainingAllRegions.txt
    fi
done
