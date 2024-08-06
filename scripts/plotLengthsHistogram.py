import numpy as np
import matplotlib.pyplot as plt
file_obj = open("/data/wrayva/output/extract_regions/regionB2/lengths.txt", "r")

# reading the data from the file
file_data = file_obj.read()

# splitting the file data into lines
lines = file_data.splitlines()
x=[]
for i in range(0, len(lines)):
    print(lines[i].split('\t'))
    print(int(lines[i].split('\t')[3]))
    value=int(lines[i].split('\t')[3])
    if value>0:
        x.append(value)
    #x.append(int(lines[i].split('\t')[3]))

fig, ax = plt.subplots(1, 1)
ax.hist(x)

# Set title
ax.set_title("Distribution of Lengths of Region B")

# adding labels
ax.set_xlabel('Length of Region B')
ax.set_ylabel('Frequency')

plt.show()
plt.savefig('/data/wrayva/output/extract_regions/regionB2/plots/regionB2lengthsHist.png')
