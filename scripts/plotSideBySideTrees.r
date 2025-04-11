library('TreeDist')
library(ape)

x=1
fileContainingListOfTreeFileName <- paste(paste("treefilenames_excludeGap", x, sep=""), ".csv", sep="")
fileNameDF <- read.csv(fileContainingListOfTreeFileName)
treeList <- vector("list", nrow(fileNameDF))

for(i in 1:nrow(fileNameDF)) {
    treeList[[i]] <- ape::read.tree(paste("../treefilecopies/", fileNameDF[i, ], sep=""))
}



jpeg('/data/Phillippy2/projects/acro_comparisons/hprc/distal/extract_regions/DJ/rfdist/plots/rfdistplot_sideBySide.jpg')
cophyloplot(treeList[1], treeList[2], assoc = NULL, use.edge.length = FALSE, space = 0,
       length.line = 1, gap = 2, type = "phylogram", rotate = FALSE,
       col = par("fg"), lwd = par("lwd"), lty = par("lty"),
       show.tip.label = TRUE, font = 3)

cophyloplot(treeList[1], treeList[2])

dev.off()