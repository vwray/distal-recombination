scriptDir=/data/wrayva/scripts


sbatch --time=5:00:00 --mem=32g --cpus-per-task=5 ${scriptDir}/run_minimap.sh /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa $regionDir/chm13chr22${region}.fa $regionDir/paf/bonobo_minimap_${region}_aln.paf

cd /data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/DJ/rfdist4

module load minimap2


#extract a sample from window 1-5000
module load samtools

samtools faidx dj_aln_1-5000.aln distal_HG00099_chr13_haplotype1-0000019:2235491-2579860:1-5000 > minimap/sample1.fna

minimap2 -cx asm20 /data/Phillippy2/projects/acro_comparisons/refs/chm13/distal_bits/chr22.distal.upper.fa minimap/sample1.fna > minimap/sample1.paf

#grab the 8th value from the paf file. this is the equivalent position in chm13 chr 22