region=DJ
extractDir=/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions
regionDir=${extractDir}/${region}
scriptDir=/data/wrayva/gitRepos/distal-recombination/scripts
#/data/wrayva/scripts
genomeDir=/data/Phillippy2/projects/hprc-assemblies/assemblies-v3
#sequencesDir=/data/wrayva/output/sequences
#subsetDir=$regionDir/50random
msaDir=${regionDir}/rfdist4
testDir=$msaDir/test

testMsaFile=$testDir/testmsa.fna

module load iqtree
cd $testDir
iqtree2 -s $testMsaFile



"sequences contain more than 50% gaps/ambiguity"

cat ${testMsaFile}.log | grep "identical sequences (see below)" #"sequences contain more than 50% gaps/ambiguity"

cp $testMsaFile copy_msa_file.fna
rm testmsa.fna*
mv copy_msa_file.fna $testMsaFile



cd $msaDir
for file in `ls *.aln.log`; do
    range=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
    newRange=$((`echo $range | cut -d ":" -f1` + 4448073))
    numberOfGapSeqs=`cat $file | grep "sequences contain more than 50% gaps/ambiguity" | cut -d " " -f2`
    if [ $numberOfGapSeqs -ge 1 ]; then
        #echo "region $range has $numberOfGapSeqs sequences containing more than 50% gaps/ambiguity" >> problemSequencesOrRegions.txt
        echo ${file%.log} >> excludeList.txt
    fi
    numberOfIdentSeqs=`cat $file | grep "identical sequences (see below) will be ignored for subsequent analysis" | cut -d " " -f2 | head -n1`
    if [ $numberOfIdentSeqs -ge 1 ]; then
        echo "$file has $numberOfIdentSeqs identical sequences" #>> problemSequencesOrRegions.txt
    fi
done

cd $msaDir/treefilecopies_excludeGap
cd $msaDir
#echo "TreeFileNames" >> $msaDir/treefilenames_excludeGap.csv
#echo "Positions" >> $msaDir/treefilepositions_excludeGap.csv
cp $msaDir/treefilenames.csv $msaDir/treefilenames_excludeGap.csv
cp $msaDir/treefilepositions.csv $msaDir/treefilepositions_excludeGap.csv
for file in `cat excludeList.txt`; do
    range=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
    second=`echo $range | cut -d '-' -f2`
    cat treefilenames_excludeGap.csv | grep $range
    cat treefilepositions_excludeGap.csv | grep $range
    sed --in-place "/$range/d" $msaDir/treefilenames_excludeGap.csv
    sed --in-place "/$second/d" $msaDir/treefilepositions_excludeGap.csv
done

for file in `cat $msaDir/treefilenames.csv`; do
    range=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
    second=`echo $range | cut -d '-' -f2`
    isInFile=$(cat excludeList.txt | grep -c $range)
    if [ $isInFile -eq 0 ]; then
        echo $file >> $msaDir/treefilenames_excludeGap.csv
        echo $second >> $msaDir/treefilepositions_excludeGap.csv
    else
        echo "$range will be excluded for gaps"
    fi
done

#put things in separate files when cross over gap
a=1
for file in `cat $msaDir/treefilenames.csv`; do
    range=`echo $file | cut -d '_' -f3 | cut -d '.' -f1`
    second=`echo $range | cut -d '-' -f2`
    isInFile=$(cat excludeList.txt | grep -c $range)
    if [ $isInFile -eq 0 ]; then
        echo $file >> "$msaDir/excludeGap/treefilenames_excludeGap${a}.csv"
        echo $second >> "$msaDir/excludeGap/treefilepositions_excludeGap${a}.csv"
    else
        echo "$range will be excluded for gaps"
        a=$((a+1))
    fi
done

