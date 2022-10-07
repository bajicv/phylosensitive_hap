## this script:
# calculates pair-wise distances from FASTA and visualizes them

## this script requires: 
# the following R packages:
library(tidyverse)
library(ape)

# loading the following data:
ape_fasta <- read.FASTA("KGP_fasta_all_subset.fasta", type = "DNA") #FASTA subset (output from creating_FASTA_subset_xyz.R)
metainfo <- read.table("metainfo_OG.txt", header=T, sep=" ")
rhb_partition <- read.csv("rhb_african_level10_singletons.csv", header = T, sep = ";") #rhierBAPS partition for merging with pariwise distances results

#################################################

# calculate pairwise distances
res <- dist.dna(ape_fasta, model = "raw", variance = FALSE,
                gamma = FALSE, pairwise.deletion = FALSE,
                base.freq = NULL, as.matrix = TRUE) 

# convert DNAbin into long format table for easier handling
res_table <- res %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  gather("Var2", "value", -Var1)

# prepare df containing info for haplogroups & rhierBAPS groups
# fuse/filter with rhierBAPS results
metainfo_partition <- left_join(rhb_partition, metainfo, by = c("Isolate" = "SampleID"))

# create Yamamoto haplogroups >> substring first letters of haplogroups and add this vector into newly created column of metainfo
metainfo_partition$SC_haplogroups <- substr(metainfo_partition$Haplogroup, 1, 1)

# create Rishishwar haplogroups >> as seen above except for haplogroups "L", where first 2 letters are chosen
for (i in 1:length(metainfo_partition$Haplogroup)) {
  if (startsWith(metainfo_partition$Haplogroup[i], "L") == T) {
    metainfo_partition$SCL_haplogroups[i] <- substr(metainfo_partition$Haplogroup[i], 1, 2)
  } else {
    metainfo_partition$SCL_haplogroups[i] <- substr(metainfo_partition$Haplogroup[i], 1, 1)
  }
}

# add to the results table of pairwise distances the info about which groupings each part of the pair is in
# Yamamoto haplogroups (so haplogroup grouping for each pair is known)
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var1" = "Isolate")) #combine for first variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'SC_haplogroups'] <- 'Var1_SC' #rename
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var2" = "Isolate")) #combine for second variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'SC_haplogroups'] <- 'Var2_SC'#rename

# Rishishwar haplogroups (so haplogroup grouping for each pair is known)
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var1" = "Isolate")) #combine for first variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'SCL_haplogroups'] <- 'Var1_SCL' #rename
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var2" = "Isolate")) #combine for second variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'SCL_haplogroups'] <- 'Var2_SCL'#rename

# rhierBAPS level 2 group (so rhierBAPS grouping for each pair is known)
res_table_groupings <- left_join(res_table, metainfo_partition, by = c("Var1" = "Isolate")) #combine for first variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'level.2'] <- 'Var1_level.2' #rename
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var2" = "Isolate")) #combine for second variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'level.2'] <- 'Var2_level.2'#rename

# rhierBAPS level 3 group (so rhierBAPS  grouping for each pair is known)
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var1" = "Isolate")) #combine for first variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'level.3'] <- 'Var1_level.3' #rename
res_table_groupings <- left_join(res_table_groupings, metainfo_partition, by = c("Var2" = "Isolate")) #combine for second variable/part of pair
names(res_table_groupings)[names(res_table_groupings) == 'level.3'] <- 'Var2_level.3'#rename

# calculate mean pairwise distance within groups 

# yamamoto
SC_mean_pairwise <- res_table_groupings %>% #filter all pairs wich were assigned to same group
  filter(Var1_SC == Var2_SC) %>%
  group_by(Var1_SC) %>% #group by all unique groups
  summarise_at(vars(value), list(name = mean)) #calculate means

# rishishwar
SCL_mean_pairwise <- res_table_groupings %>% #filter all pairs wich were assigned to same group
  filter(Var1_SCL == Var2_SCL) %>%
  group_by(Var1_SCL) %>% #group by all unique groups
  summarise_at(vars(value), list(name = mean)) #calculate means

# rhb lev 2
rhb_level2_mean_pairwise <- res_table_groupings %>%
  filter(Var1_level.2 == Var2_level.2) %>% #filter all pairs wich were assigned to same group
  group_by(Var1_level.2) %>% #group by all unique groups
  summarise_at(vars(value), list(name = mean)) #calculate means

# rhb lev 3
rhb_level3_mean_pairwise <- res_table_groupings %>% #filter all pairs wich were assigned to same group
  filter(Var1_level.3 == Var2_level.3) %>%
  group_by(Var1_level.3) %>% #group by all unique groups
  summarise_at(vars(value), list(name = mean)) #calculate means

#subset per grouping method
rhb_SC_subset <- res_table_groupings %>%
  subset(Var1_SC == Var2_SC)

rhb_SCL_subset <- res_table_groupings %>%
  subset(Var1_SCL == Var2_SCL)

rhb_level2_subset <- res_table_groupings %>%
  subset(Var1_level.2 == Var2_level.2)

rhb_level3_subset <- res_table_groupings %>%
  subset(Var1_level.3 == Var2_level.3)

# color palette matching for_tree_visualization color coding
colors_SC <- c("A" = "#969696", "B" = "#f0f0f0", "C" = "#525252", 
                      "D" = "#737373", "H" = "#ffffff", "J" = "#bdbdbd", 
                      "M" = "#252525", "U" = "#d9d9d9", "L" = "#3182bd")

colors_SCL <- c("A" = "#969696", "B" = "#f0f0f0", "C" = "#525252", 
                      "D" = "#737373", "H" = "#ffffff", "J" = "#bdbdbd", "M" = "#252525", "U" = "#d9d9d9", 
                      "L3" = "#de2d26", "L4" = "#fee0d2", "L2" = "#3182bd", "L5" = "#bf812d",
                      "L1" = "#8c510a", "L0" = "#31a354")

colors_rhb2 <- c("1" = "#de2d26", "2" = "#6a51a3", "3" = "#cbc9e2", "6" = "#3182bd", 
                      "4" = "#fd8d3c", "5" = "#31a354",  "7" = "#8c510a")

colors_rhb3 <- c("1" = "#de2d26", "2" = "#fc9272", "4" = "#6a51a3", "3" = "#9e9ac8", "13" = "#cbc9e2", 
                      "12" = "#3182bd", "11" = "#9ecae1", "10" = "#deebf7", 
                      "6" = "#fd8d3c", "5" = "#fdbe85", "9" = "#31a354", "8" = "#a1d99b", 
                      "7" = "#e5f5e0", "14" = "#8c510a", "15" = "#dfc27d", "16" = "#bf812d")

#plot as violinplots
#change order of factors for matching across groupings
rhb_level2_subset$Var1_level.2 <- factor(rhb_level2_subset$Var1_level.2, levels = c("1", "2", "3", "6", "4", "5", "7")) 

rhb2_vio <- ggplot(rhb_level2_subset, aes(x= factor(Var1_level.2), fill = factor(Var1_level.2), y=value)) + 
  geom_violin() +
  scale_fill_manual(values = colors_rhb2) +
  ylim(0, 0.008) +
  geom_point(stat = "summary", fun = "mean") +
  labs(fill = "rhB level 2") +
  ylab("pairwise distance") + xlab(element_blank()) +
  theme(text = element_text(size = 25)) 
rhb2_vio 

#change order of factors for matching across groupings
rhb_level3_subset$Var1_level.3 <- factor(rhb_level3_subset$Var1_level.3, levels = c("1", "2", "4", "3", "13", "12", "11", "10", "6", "5", "9", "8", "7", "14", "15", "16")) 

rhb3_vio <- ggplot(rhb_level3_subset, aes(x= factor(Var1_level.3), fill = factor(Var1_level.3), y=value)) + 
  geom_violin() +
  scale_fill_manual(values = colors_rhb3) +
  ylim(0, 0.008) +
  geom_point(stat = "summary", fun = "mean") +
  labs(fill = "rhB level 3") +
  ylab("pairwise distance") + xlab(element_blank()) +
  theme(text = element_text(size = 25)) 
rhb3_vio

#plot
#change order of factors for matching across groupings
rhb_yama_subset$Var1_SC <- factor(rhb_yama_subset$Var1_SC, levels = c("A", "B", "C", "D", "H", "J", "M", "U", "L")) 

SC_vio <- ggplot(rhb_SC_subset, aes(x= factor(Var1_SC), fill = factor(Var1_SC), y=value)) + 
  geom_violin() +
  scale_fill_manual(values = colors_SC) +
  ylim(0, 0.008) +
  geom_point(stat = "summary", fun = "mean") +
  labs(fill = "SC+L") +
  ylab("pairwise distance") + xlab(element_blank()) +
  theme(text = element_text(size = 25)) 
SC_vio

#change order of factors for matching across groupings
rhb_SCL_subset$Var1_SCL <- factor(rhb_rish_subset$Var1_SCL, levels = c("A", "B", "C", "D", "H", "J", "M", "U", "L3", "L4", "L2", "L5", "L1", "L0")) 

SCL_vio <- ggplot(rhb_SCL_subset, aes(x= factor(Var1_SCL), fill = factor(Var1_SCL), y=value)) + 
  geom_violin() +
  scale_fill_manual(values = colors_SCL) +
  ylim(0, 0.008) +
  geom_point(stat = "summary", fun = "mean") +
  labs(fill = "SC+L") +
  ylab("pairwise distance") + xlab(element_blank()) +
  theme(text = element_text(size = 25)) 
SCL_vio

