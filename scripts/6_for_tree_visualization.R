## this script:
# subsets a FASTA file based on the superpopulation "AFR" from the 1000 Genomes Project (1KGP)
# creates a FASTA file used as input for various downstream analyses

## this script requires: 
# the following R packages:
library(tidyverse)
library(ggtree)
library(phytools)

# loading the following data:
iqtree_subset <- read.newick("KGP_fasta_visualsubset_RSRS_iqtree_tree.treefile") #output iqtree
subset_table <- read.csv("visualizing_tree.csv", header = T, sep = ";") #table containing info about the chosen representative haplogroups + grouping info
colors_metainfo <- read.table("Baja_colors_for_plot.txt")

#################################################

# prepare informative tip labels to be haplogroups instead of sample IDs
labels <- as.data.frame(iqtree_subset$tip.label) # extract order of IDs for tip labels
labels <- left_join(labels, subset_table, by = c("iqtree_subset$tip.label" = "Sample")) # order haplogroup labels accordingly

#set haplogroup labels as new tip labels
iqtree_subset$tip.label <- labels$Group.lables

# create tree with haplogroup as tip labels
gg_subset <- ggtree(iqtree_subset, layout = "rectangular", branch.length = 'none')

gg <- gg_subset +
        geom_tiplab()
gg

# use handmade list of colors and matching levels across different groupings
colors_metainfo <- colors_metainfo + 100
colors_metainfo <- mutate_if(colors_metainfo, is.numeric, as.factor)

# choose color palette
my_26_colors <- c(
  "#3182bd",
  "#9ecae1",
  "#deebf7",
  "#31a354",
  "#a1d99b",
  "#e5f5e0",
  "#fd8d3c",
  "#fdbe85",
  "#6a51a3",
  "#9e9ac8",
  "#cbc9e2",
  "#de2d26",
  "#fc9272",
  "#fee0d2",
  "#ffffff",
  "#f0f0f0",
  "#d9d9d9",
  "#bdbdbd",
  "#969696",
  "#737373",
  "#525252",
  "#252525",
  "#8c510a",
  "#bf812d",
  "#dfc27d",
  "#ffffff"
)

# plot tree + columns/boxes containing grouping correspondence across groupings (encoded by colors)
gg <- gheatmap(gg, colors_metainfo, width = .4, offset = 2.5, 
              colnames = T) +
  scale_fill_manual(values = my_26_colors) +
  theme(legend.position = "none")
gg

# save plot
ggsave("tree_with_colors.pdf", plot = last_plot())
