## this script:
# performs correspondence analysis
# and visualized it

## this script requires: 
# the following R packages:
library("tidyverse") #data handling
library("FactoMineR") #CA
library("factoextra") #visualization

# loading the following data:
metainfo <- read.table("metainfo_OG.txt", header=T, sep=" ")
rhb_partition <- read.csv("rhb_all_level7_singletons.csv", header = T, sep = ",") #rhierBAPS partition for merging with pariwise distances results

#################################################

# fuse/filter with rhierBAPS results
metainfo_partition <- left_join(rhb_partition, metainfo, by = c("Isolate" = "SampleID"))

# create Yamamoto haplogroups by extracting first letter of each haplogroup 
metainfo_partition$SC_haplogroups <- substr(metainfo_partition$Haplogroup, 1, 1)

# create Rishishwar haplogroups by extracting first 2 letters for L haplogroups and only first letter for other haplogroup 
for (i in 1:length(metainfo_partition$Haplogroup)) {
  if (startsWith(metainfo_partition$Haplogroup[i], "L") == T) {
    metainfo_partition$SCL_haplogroups[i] <- substr(metainfo_partition$Haplogroup[i], 1, 2)
  } else {
    metainfo_partition$SCL_haplogroups[i] <- substr(metainfo_partition$Haplogroup[i], 1, 1)
  }
}


# create frequency tables for Yamamoto, Rishishwar and rhierBAPS
raw_ftable <- as.data.frame.matrix(table(metainfo_partition$Population.code, metainfo_partition$Haplogroup)) #define object w/table parameters for simple calling
SC_ftable <- as.data.frame.matrix(table(metainfo_partition$Population.code, metainfo_partition$SC_haplogroups)) #define object w/table parameters for simple calling
SCL_ftable <- as.data.frame.matrix(table(metainfo_partition$Population.code, metainfo_partition$SCL_haplogroups)) #define object w/table parameters for simple calling
rhb_level2_ftable <- as.data.frame.matrix(table(metainfo_partition$Population.code, metainfo_partition$level.2)) #define object w/table parameters for simple calling
rhb_level3_ftable <- as.data.frame.matrix(table(metainfo_partition$Population.code, metainfo_partition$level.3)) #define object w/table parameters for simple calling
rhb_level6_ftable <- as.data.frame.matrix(table(metainfo_partition$Population.code, metainfo_partition$level.6)) #define object w/table parameters for simple calling

# CA
options(ggrepel.max.overlaps = Inf) #increase accepted overlap for plotting

ca_raw <- CA(raw_ftable, graph = FALSE)
ca_SC <- CA(SC_ftable, graph = FALSE)
ca_SCL <- CA(SCL_ftable, graph = FALSE)
ca_rhb_level2 <- CA(rhb_level2_ftable, graph = FALSE)
ca_rhb_level3 <- CA(rhb_level3_ftable, graph = FALSE)
ca_rhb_level6 <- CA(rhb_level6_ftable, graph = FALSE)

# plot 
colforCAplot <- c("#AC88FF", "#AC88FF", "#00A5FF", "#00B81F", "#BB9D00", "#00B81F", "#00B81F",
                  "#F8766D", "#AC88FF", "#BB9D00", "#BB9D00", "#00A5FF", "#AC88FF", "#BB9D00",
                  "#00A5FF", "#00B81F", "#00B81F", "#AC88FF", "#AC88FF", "#F8766D", "#F8766D",
                  "#00A5FF", "#F8766D", "#00A5FF", "#BB9D00", "#AC88FF")

# CA plot
plot <- fviz_ca_row(ca_raw, col.col = "black", col.row = as.list(colforCAplot), 
                    labelsize = 5, pointsize = 5, repel = T)
plot

# plot with arrows
finalplot <- fviz_ca_biplot(ca_rhb_level6, map ="colgreen", col.row = as.list(colforCAplot), col.col = "grey", 
                            labelsize = 5, pointsize = 5,
                            arrow = c(FALSE, TRUE), repel = TRUE) 
finalplot

ggsave("CA_6.pdf", finalplot, limitsize = FALSE)

