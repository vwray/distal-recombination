#!/usr/bin/env Rscript
# args[1] is the tree file
# args[2] is the output file for color by chromosome
# args[3] is the output file for color by population
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ggtree)
tree <- read.tree(args[1])
chrs <- list(chr13=c(tree$tip.label[grep("chr13", tree$tip.label)]),
             chr14=c(tree$tip.label[grep("chr14", tree$tip.label)]),
             chr15=c(tree$tip.label[grep("chr15", tree$tip.label)]),
             chr21=c(tree$tip.label[grep("chr21", tree$tip.label)]),
             chr22=c(tree$tip.label[grep("chr22", tree$tip.label)]),
             chr4=c(tree$tip.label[grep("chr4", tree$tip.label)]),
             bonobo=c(tree$tip.label[grep("mPanPan", tree$tip.label)])
)
grouped_tree1 <- groupOTU(tree, chrs,
                       group_name = "chromosome")
ggtree(grouped_tree1, aes(color=chromosome), size=3, layout='circular') +
  theme(legend.position=c(1,.98), legend.key.size = unit(1, 'cm'),plot.margin=margin(150, 150, 150, 150), legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(2))) +
  geom_tiplab(align=TRUE, linesize=.5, size=6, fontface='bold')

ggsave(args[2], width = 60, height = 60, units = "cm", limitsize = FALSE)

pops <- list(EUR=c(tree$tip.label[grep("EUR", tree$tip.label)]),
             EAS=c(tree$tip.label[grep("EAS", tree$tip.label)]),
             AMR=c(tree$tip.label[grep("AMR", tree$tip.label)]),
             AFR=c(tree$tip.label[grep("AFR", tree$tip.label)]),
             SAS=c(tree$tip.label[grep("SAS", tree$tip.label)])
)
grouped_tree2 <- groupOTU(tree, pops,
                       group_name = "population")
ggtree(grouped_tree2, aes(color=population), size=3, layout='circular') +
  theme(legend.position=c(1, .98), legend.key.size = unit(1, 'cm'),plot.margin=margin(150, 150, 150, 150), legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(2))) +
  geom_tiplab(align=TRUE, linesize=.5, size=6, fontface='bold')

ggsave(args[3], width = 60, height = 60, units = "cm", limitsize = FALSE)
