library('TreeDist')
#if the above does not work, then install it
#install.packages('TreeDist')
#filename <- "/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/DJ/50random/mafft/mafft_DJ_1-10000_copy.aln.treefile"
#tree1 <- ape::read.tree(filename)
#filename2 <- "/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/DJ/50random/mafft/mafft_DJ_10000-20000_copy.aln.treefile"
#tree2 <- ape::read.tree(filename2)
#distance <- TreeDistance(tree1, tree2)

#test
#tree1 <- ape::read.tree(text = '(A, ((B, (C, (D, E))), ((F, G))));')
#tree2 has two extra leaves, H and I
#tree2 <- ape::read.tree(text = '(A, ((B, (C, (D, (H, I)))), ((F, G), E)));')
#distance <- TreeDistance(tree1, tree2)
#0.3417448
#tree3 is same as tree2, but with leaves H and I removed
#tree3 <- ape::read.tree(text = '(A, ((B, (C, D)), ((F, G), E)));')
#distance2 <- TreeDistance(tree1, tree3)
#also 0.3417448


#report split info, clustering info
#ClusteringInfoDistance(tree1, tree2, reportMatching = TRUE)

fileContainingListOfTreeFileName <- "treefilenames.csv" #"treefilenames_odd.csv"
fileNameDF <- read.csv(fileContainingListOfTreeFileName)
treeList <- vector("list", nrow(fileNameDF))

for(i in 1:nrow(fileNameDF)) {
    treeList[[i]] <- ape::read.tree(paste("treefilecopies/", fileNameDF[i, ], sep=""))
}

distanceList <- vector("list", (length(treeList) - 1))
for(i in 1:(length(treeList) - 1)){
    distanceList[[i]] <- TreeDistance(treeList[[i]], treeList[[i+1]])
}

coordsfile <- "treefilepositions.csv" #"treefilepositions_odd.csv"
coordsDF <- read.csv(coordsfile)
coords <- as.list(coordsDF)
#c(10000,20000,30000,40000)
coords <- coordsDF$Positions

newPositions <- coords[1:length(coords)-1] + 4448073

jpeg('/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/regionA/chunk/plots/rfdistplot.jpg')
plot(coords[1:length(coords)-1], distanceList, type="l", col="red", main="RF Distance Between Adjacent Trees: Region A",
   xlab="Position", ylab="RF Distance", pch=19)
#plot(newPositions, distanceList, type="l", col="red", main="RF Distance Between Adjacent Trees (5kb chunks, 2.5kb overlap, \n Excluding Gaps)",
#   xlab="Position on CHM13 Chromosome 22", ylab="RF Distance", pch=19)
dev.off()
