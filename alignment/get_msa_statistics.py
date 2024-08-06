import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as p

fastaObject2 = pysam.FastaFile('/data/wrayva/output/extract_regions/DJ/50random/mafft/mafft_DJ.aln')

length_of_msa = len(fastaObject2.fetch(reference=fastaObject2.references[1]))
rows=5
cols=length_of_msa
counts= [[0 for i in range(cols)] for j in range(rows)]

for i in range(len(fastaObject2.references)):
    print(fastaObject2.references[i])
    sequence = fastaObject2.fetch(reference=fastaObject2.references[i])
    for j in range(length_of_msa):
        if(sequence[j]=='A' or sequence[j]=='a'):
            counts[0][j]+=1
        if(sequence[j]=='C' or sequence[j]=='c'):
            counts[1][j]+=1
        if(sequence[j]=='T' or sequence[j]=='t'):
            counts[2][j]+=1
        if(sequence[j]=='G' or sequence[j]=='g'):
            counts[3][j]+=1
        if(sequence[j]=='-' or sequence[j]=='N'):
            counts[4][j]+=1

x=[]
mut_counts=[]
for j in range(length_of_msa):
    #take max of counts[length_of_msa] as consensus, get count of ones not in consensus
    index = p.Series([counts[0][j], counts[1][j], counts[2][j], counts[3][j], counts[4][j]]).idxmax()
    # if max id is 4 indicating '-' then continue
    if(index==4):
        continue
    numMutations = 0
    for i in range(4):
        if(i!=index):
            numMutations += counts[i][j]
    mut_counts.append(numMutations)
    x.append(j)



divisor = float(len(fastaObject2.references))
mut_percents = [float(x)/divisor for x in mut_counts]

ave_percent = np.average(mut_percents)
ave_count = np.average(mut_counts)
print(f"total sequences: {divisor} ave mutation count: {ave_count} ave mutation percent: {ave_percent}")

plt.clf()

plt.scatter(x, mut_percents, s=1, linewidth=0 )
#plt.plot(x, mut_percents)
plt.xlabel("Position along alignment")
plt.ylabel("Proportion of sequences with mutation")
plt.title("Proportion of sequences with mutation at each position in DJ")
plt.savefig('/data/wrayva/output/plots/dj_msa_avg_mut_skip_gap_50_random.pdf')


#In the DJ, the average number of sequences with mutations at each position is 165 out of a total 267 sequences, which is ~61.8% of sequences. In region A, the average number of sequences with mutations at each position is less than 1 (0.2) out of a total 291 sequences, which is ~0.0007% of sequences. I have plotted the percentage of sequences with mutations at each position. Note that gaps are not being counted as mutations (should this be handled differently?) so this is may be skewing the numbers lower.

#region A ignoring positions with majority gap: total sequences: 272.0 ave mutation count: 0.40228945216680295 ave mutation percent: 0.0014790053388485402
#dj ignoring positions with majority gap: total sequences: 267.0 ave mutation count: 1.6137193645047412 ave mutation percent: 0.006043892750954089
#region B ignoring positions with majority gap: total sequences: 326.0 ave mutation count: 1.1725612090496418 ave mutation percent: 0.0035968135246921534
#region C: total sequences: total sequences: 233.0 ave mutation count: 0.6451431425897096 ave mutation percent: 0.002768854689226221

#testing
m = max(counts[0][j], counts[1][j], counts[2][j], counts[3][j])
print(m)
index = p.Series([counts[0][j], counts[1][j], counts[2][j], counts[3][j]]).idxmax()
print(index)
numMutations = 0
for i in range(4):
    if(i!=index):
        numMutations += counts[i][j]





sequence = fastaObject.fetch(reference=fastaObject2.references[1])
for j in range(length_of_msa):
    if(sequence[j]=='A' or sequence[j]=='a'):
        counts[0][j]+=1

sequence = fastaObject1.fetch(reference='distal_NA20905_chr22_haplotype2-0000055:2701068-3046855')



sequence = "ACTG-AAA"
length_of_msa=len(sequence)
rows=5
cols=length_of_msa
counts= [[0 for i in range(cols)] for j in range(rows)]
for j in range(length_of_msa):
    if(sequence[j]=='A' or sequence[j]=='a'):
        counts[0][j]+=1
    if(sequence[j]=='C' or sequence[j]=='c'):
        counts[1][j]+=1
    if(sequence[j]=='T' or sequence[j]=='t'):
        counts[2][j]+=1
    if(sequence[j]=='G' or sequence[j]=='g'):
        counts[3][j]+=1
    if(sequence[j]=='-' or sequence[j]=='N'):
        counts[4][j]+=1



'''

[E::fai_retrieve] Failed to retrieve block: unexpected end of file
Traceback (most recent call last):
  File "<stdin>", line 9, in <module>
  File "pysam/libcfaidx.pyx", line 317, in pysam.libcfaidx.FastaFile.fetch
ValueError: failure when retrieving sequence on 'distal_NA20905_chr22_haplotype2-0000055:2701068-3046855'
'''
