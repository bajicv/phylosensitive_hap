## this script:
# visualizes the partition dataframe (output from rhierBAPS)
# using a phylogenetic tree (output from running iqtree) 
# IMPORTANT: the annotation only works if the sample IDs of the partition and tree match, 
# as ensured in the script creating each of the input FASTAs 

## this script requires: 
# the following R packages:
library(tidyverse)
library(phytools)
library(ggtree)

# loading the following data:
rhb_partition <- read.csv("rhb_african_level10_singletons.csv", header = T, sep = ";") #partition dataframe (from rhierBAPS results (output from rhierBAPS.R)
iqtree <- phytools::read.newick("KGP_fasta_RSRS_african_iqtree_tree.treefile") #treefile (output from running iqtree in commandline)
metainfo <- read.table("metainfo_OG.txt", header=T, sep=" ") #metainfo file (output from metainfo.R) for labelling tree with metainfo

#################################################

# choose tree layout
gg <- ggtree(iqtree, layout = "circular", branch.length = 'none')

metainfo_partition <- left_join(rhb_partition, metainfo, by = c("Isolate" = "SampleID")) %>%
  unite("meta", c(Population.name, Haplogroup))

# plot the annotated tree
gg <- gg %<+% metainfo_partition #add df containing new tip labels to tree
gg <- gg + geom_tippoint(aes(color = factor(level.3))) + 
  geom_tiplab(aes(label = meta), size = 1, offset = 1) +
  theme(legend.position="top") +
  theme(plot.margin = unit(c(15, 15, 15, 15), "mm")) 
gg

# save tree
ggsave("AFR_tree_level3.pdf", gg, width = 50, height = 50, units = "cm", limitsize = FALSE)
