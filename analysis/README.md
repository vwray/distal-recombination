This document contains details on the analysis I have run as part of this project.

## Align distal bits with CHM13 with minimap:
```sh
minimap2 -x asm20 /data/Phillippy/references/T2T-CHM13/chm13v2.0.fa /data/Phillippy2/projects/acro_comparisons/hprc/distal/distalbits.fasta > /data/Phillippy2/projects/acro_comparisons/hprc/chm13_to_distalbits.paf
```

## Create histogram of distal bit lengths:
This bash code first creates a file `distalbits_lineLengths.txt` containing the line lengths of each line in the fasta file. Then we remove the line lengths corresponding to sequence identifier headings to create `distalbits_modifiedlineLengths.txt`.
```bash
awk '{ print length }' distalbits.fasta >> distalbits_lineLengths.txt

i=1
for line in `cat distalbits_lineLengths.txt`; do
    if [[ $(($i % 2)) -eq 0 ]]; then
        echo $line >> distalbits_modifiedlineLengths.txt
    fi
    ((i++))
done
```
This python code creates a histogram from the line lengths, which are the lengths of the distal regions:
```python
import numpy as np
import matplotlib.pyplot as plt
file_obj = open("/path/distalbits_modifiedlineLengths.txt", "r")

# create a list containing each length from the file, which are all on separate lines in the file
file_data = file_obj.read()
lines = file_data.splitlines()
x=[]
for i in range(0, len(lines)):
    x.append(int(lines[i]))

fig, ax = plt.subplots(1, 1)
ax.hist(x)

# Add title and axis labels
ax.set_title("Distribution of Lengths of Distal Bits")=
ax.set_xlabel('Length of distal bit')
ax.set_ylabel('Frequency')

plt.show()
plt.savefig('/path/distalbits_lengthsHistogram.png')
```

## Create %GC Content Plot:
This Python script creates a plot of the %GC content of a random selection of 40 of the sequences.
```python
from Bio.SeqUtils import gc_fraction
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random

def compute_gc_content(record, window_size, step_size, ax, offset):
    """Plot windowed GC content on a designated Matplotlib ax."""
    #Computes GC content based on window and step size
    #Add the offset in order to stack sequence plots
    y = [100.0*gc_fraction(record[i : i + window_size])+offset for i in range(0,len(record) - window_size, step_size)]
    x = [i for i in range(0,len(record) - window_size, step_size)]
    y2 = [offset for i in range(0,len(record) - window_size, step_size)]
    ax.fill_between(x, y, y2, alpha=0.5)

#Opens the fasta file
fastaObject = pysam.FastaFile('/data/wrayva/output/wfmash_on_trimmed_to_chm13/query.fa')

#Set up the title and axis labels
fig, ax = plt.subplots(1, 1)
ax.set_ylim(bottom=0,top=4000)
ax.set_ylabel("GC(%)")
ax.set_xlabel("Distal Region Position")
plt.yticks([])
plt.title("%GC Content for 40 Random Distal Regions\n(Window Size 1000, Step Size 100)")

#Chooses 40 random sequences in the fasta file
random_indices = set()
while(len(random_indices)< 40):
    random_indices.add(random.randint(0,len(fastaObject.references)-1))
random_indices_list = list(random_indices)

max_length = 0
for i in range(len(fastaObject.references)):
    j=random_indices_list[i]
    region = fastaObject.fetch(region=fastaObject.references[j])
    #Add the sequence name label to the plot
    ax.annotate(text=fastaObject.references[j][7:24]+fastaObject.references[j][30], xy=(-1400000,i*100+10), fontsize=0.5)
    compute_gc_content(region, 1000, 100, ax, i*100)
    if len(region) > max_length:
        max_length = len(region)

ax.set_xlim([-1500000,max_length+100000])
fig.savefig('/data/wrayva/output/plots/gc_plot_1.pdf')
```
