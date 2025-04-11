library('TreeDist')

jpeg('/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/DJ/rfdist/plots/rfdistplot_skipGapsScatter.jpg')
par(xpd=TRUE)
for (x in c(1,3,9,12)) {
    #print(paste(paste("treefilenames_excludeGap", x, sep=""), ".csv", sep=""))



    fileContainingListOfTreeFileName <- paste(paste("treefilenames_excludeGap", x, sep=""), ".csv", sep="")
    fileNameDF <- read.csv(fileContainingListOfTreeFileName)
    treeList <- vector("list", nrow(fileNameDF))

    for(i in 1:nrow(fileNameDF)) {
        treeList[[i]] <- ape::read.tree(paste("../treefilecopies/", fileNameDF[i, ], sep=""))
    }

    distanceList <- vector("list", (length(treeList) - 1))
    for(i in 1:(length(treeList) - 1)){
        distanceList[[i]] <- TreeDistance(treeList[[i]], treeList[[i+1]])
    }

    coordsfile <- paste(paste("treefilepositions_excludeGap", x, sep=""), ".csv", sep="")
    coordsDF <- read.csv(coordsfile)
    coords <- as.list(coordsDF)
    #c(10000,20000,30000,40000)
    coords <- coordsDF$Positions

    newPositions <- coords[1:length(coords)-1] + 4448073

    #plot(coords[1:length(coords)-1], distanceList, type="l", col="red", main="RF Distance Between Adjacent Trees: Region A",
    #xlab="Position", ylab="RF Distance", pch=19)
    if (x == 1) {
        plot(newPositions, distanceList, type="p", col="red", main="RF Distance Between Adjacent Trees within DJ \n (5kb chunks, 2.5kb overlap, Excluding Gaps)",
            xlab="Position on CHM13 Chromosome 22", ylab="RF Distance", pch=19, xlim = c(4448073,4813073), ylim=c(0,1))
        print("Initial coords")
    } else {
        points(newPositions, distanceList, type="p", col="red", pch=19, xpd = TRUE)
        print("Adding coords from file")
        print(x)
    }

}
dev.off()
