This document contains details on the analysis I have run as part of this project.

## Align distal bits with CHM13 with minimap:
```
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
