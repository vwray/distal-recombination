from Bio.SeqUtils import gc_fraction
import pysam
import matplotlib.pyplot as plt
import numpy as np

with open('/data/wrayva/output/plots/wiggleFile4.wig', 'w') as f:
    def add_to_wiggle_file(record, window_size, step_size, chromosome, name, start):
        f.write('track type="wiggle_0" name='+ name + '\n')
        f.write('fixedStep  chrom=chr'+ chromosome + ' start=' + start + ' step=100 span=1000\n')
        for i in range(0,len(record) - window_size, step_size):
            f.write(str(100.0*gc_fraction(record[i : i + window_size])))
            f.write('\n')


    fastaObject2 = pysam.FastaFile('/data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa')
    for i in range(len(fastaObject2.references)):
        region2 = fastaObject2.fetch(region=fastaObject2.references[i])
        chromosome = fastaObject2.references[i][18:20]
        name = fastaObject2.references[i][7:]
        print(chromosome)
        print(name)
        fullname = fastaObject2.references[i]
        t_filename = "/data/wrayva/output/seqtk_telo/output_rdna_trim_"+fullname+".txt"
        print(t_filename)
        telomere_file = open(t_filename, "r")
        start = telomere_file.readline().split('\t')[2]
        add_to_wiggle_file(region2, 1000, 100, chromosome, name, start)



fastaObject2 = pysam.FastaFile('/data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa')
region2 = fastaObject2.fetch(region=fastaObject2.references[0])
chromosome = fastaObject2.references[0][18:20]
name = fastaObject2.references[0][7:]
fullname = fastaObject2.references[0]
t_filename = "/data/wrayva/output/seqtk_telo/output_rdna_trim_"+fullname+".txt"
print(t_filename)
telomere_file = open(t_filename, "r")
start = telomere_file.readline().split('\t')[2]





from Bio.SeqUtils import gc_fraction
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random


def compute_gc_content(record, window_size, step_size, ax, offset):
    """Plot windowed GC content on a designated Matplotlib ax."""
    y = [100.0*gc_fraction(record[i : i + window_size])+offset for i in range(0,len(record) - window_size, step_size)]
    x = [i for i in range(0,len(record) - window_size, step_size)]
    y2 = [offset for i in range(0,len(record) - window_size, step_size)]
    #np.arange(len(record) - window_size) + 25
    ax.fill_between(x, y, y2, alpha=0.5)

    #ax.plot(x,y)


fastaObject = pysam.FastaFile('/data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa')
#print("number of entries in fasta file:", len(fastaObject.references))
#391 entries in query.fa

fig, ax = plt.subplots(1, 1)
ax.set_ylim(bottom=0,top=39100)
#ax.set_xlim(left=-2800000)
ax.set_ylabel("GC(%)")
ax.set_xlabel("Distal Region Position")
plt.yticks([])
plt.title("%GC Content for All Distal Regions\n(Window Size 1000, Step Size 100)")

random_indices = set()
while(len(random_indices)< 40):
    random_indices.add(random.randint(0,len(fastaObject.references)-1))

print(random_indices)
print(len(random_indices))
random_indices_list = list(random_indices)


max_length = 0

for i in range(len(fastaObject.references)):
    j=i#random_indices_list[i]
    region = fastaObject.fetch(region=fastaObject.references[j])
    ax.annotate(text=fastaObject.references[j][7:24]+fastaObject.references[j][30], xy=(-250000,i*100+10), fontsize=0.5)
    compute_gc_content(region, 1000, 100, ax, i*100)
    if len(region) > max_length:
        max_length = len(region)


ax.set_xlim([-300000,max_length+100000])

#region1 = fastaObject1.fetch(region=fastaObject1.references[0])

#fastaObject2 = pysam.FastaFile('/data/wrayva/output/trim_telo/trim_telo_rdna_trim_distal_HG04184_chr22_haplotype1-0000004.fna')
#region2 = fastaObject2.fetch(region=fastaObject2.references[0])

#fastaObject3 = pysam.FastaFile('/data/wrayva/output/trim_telo/trim_telo_rdna_trim_distal_HG02004_chr14_haplotype1-0000029.fna')
#region3 = fastaObject3.fetch(region=fastaObject3.references[0])

#compute_gc_content(region1, 1000, 100, ax, 0)
#compute_gc_content(region2, 1000, 100, ax, 100)
#compute_gc_content(region3, 1000, 100, ax, 200)

fig.savefig('/data/wrayva/output/plots/gc_18.pdf')
