## this script:
# recreates macro-haplogroups as defined by Yamamoto and Rishishwar
# recreates frequency barplots as seen in Yamamoto et al for Yamamoto, Rishishwar, and rhierBAPS groupings fro AFR subset

## this script requires: 
# the following R packages:
library(tidyverse)

# loading the following data:
metainfo <- read.table("metainfo_OG.txt", header=T, sep=" ") #metainfo file (created in metainfo.R)
rhb_partition <- read.csv("rhb_african_level3_singletons_default.csv", header = T, sep = ",") #partition dataframe (from rhierBAPS results (output from rhierBAPS.R) 

################################################# 

# subset metainfo by superpopulation == AFR to acquire samples assigned AFR
metainfo_afr <- metainfo %>% 
  subset(Superpopulation.code == "AFR") 

# create Yamamoto haplogroups by extracting first letter of each haplogroup 
metainfo_afr$SC_haplogroups <- substr(metainfo_afr$Haplogroup, 1, 1)

# create Rishishwar haplogroups by extracting first 2 letters for L haplogroups and only first letter for other haplogroup 
for (i in 1:length(metainfo_afr$Haplogroup)) {
  if (startsWith(metainfo_afr$Haplogroup[i], "L") == T) {
    metainfo_afr$SCL_haplogroups[i] <- substr(metainfo_afr$Haplogroup[i], 1, 2)
  } else {
    metainfo_afr$SCL_haplogroups[i] <- substr(metainfo_afr$Haplogroup[i], 1, 1)
  }
}

# add rhierBAPS partition to metainfo
metainfo_afr <- left_join(metainfo_afr, rhb_partition, by = c("SampleID" = "Isolate"))
metainfo_afr <- metainfo_afr %>%
  select(SampleID, level.1, level.2, level.3, Haplogroup, Population.code, Superpopulation.code, SC_haplogroups, SCL_haplogroups) #select for tidier overview
names(metainfo_afr)[names(metainfo_afr) == 'level.1'] <- 'rhierBAPS_lev1' #rename for convenience
names(metainfo_afr)[names(metainfo_afr) == 'level.2'] <- 'rhierBAPS_lev2' #rename for convenience
names(metainfo_afr)[names(metainfo_afr) == 'level.3'] <- 'rhierBAPS_lev3' #rename for convenience

# color palette matching for_tree_visualization color coding
colors_SC_tree <- c("A" = "#969696", "B" = "#f0f0f0", "C" = "#525252", 
                    "D" = "#737373", "H" = "#ffffff", "J" = "#bdbdbd", 
                    "M" = "#252525", "U" = "#d9d9d9", "L" = "#3182bd")

colors_SCL_tree <- c("A" = "#969696", "B" = "#f0f0f0", "C" = "#525252", 
                    "D" = "#737373", "H" = "#ffffff", "J" = "#bdbdbd", "M" = "#252525", "U" = "#d9d9d9", 
                    "L3" = "#de2d26", "L4" = "#fee0d2", "L2" = "#3182bd", "L5" = "#bf812d",
                    "L1" = "#8c510a", "L0" = "#31a354")

colors_rhb1_tree <- c("1" = "#de2d26", "2" ="#3182bd")

colors_rhb2_tree <- c("1" = "#de2d26", "2" = "#6a51a3", "3" = "#cbc9e2", "6" = "#3182bd", 
                      "4" = "#fd8d3c", "5" = "#31a354",  "7" = "#8c510a")

colors_rhb3_tree <- c("1" = "#de2d26", "2" = "#fc9272", "4" = "#6a51a3", "3" = "#9e9ac8", "13" = "#cbc9e2", 
                      "12" = "#3182bd", "11" = "#9ecae1", "10" = "#deebf7", 
                      "6" = "#fd8d3c", "5" = "#fdbe85", "9" = "#31a354", "8" = "#a1d99b", 
                      "7" = "#e5f5e0", "14" = "#8c510a", "15" = "#dfc27d", "16" = "#bf812d")

# raw haplogroup frquency barplot
raw_AFR_barplot <- ggplot(metainfo_afr) + 
  geom_bar(mapping = aes(x = Superpopulation.code, fill = Haplogroup, y = (..count..)/sum(..count..)), position = position_fill(reverse = TRUE)) +
  labs(fill = "raw haplogroups") +
  xlab(element_blank()) + ylab("haplogroup proportions") +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "white")
  ) +
  theme(legend.position="none") +
  coord_flip()
raw_AFR_barplot

# save raw in same way for figure layout
ggsave("frequency_barplots_raw_afr.pdf", raw_AFR_barplot, limitsize = FALSE)

# Yamamoto frequency barplot
#change order of factors for matching across groupings
metainfo_afr$SC_haplogroups <- factor(metainfo_afr$SC_haplogroups, levels = c("A", "B", "C", "D", "H", "J", "M", "U", "L")) 

#plot
SC_AFR_barplot <- ggplot(metainfo_afr) + 
  geom_bar(mapping = aes(x = Superpopulation.code, fill = SC_haplogroups, y = (..count..)/sum(..count..)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values=colors_SC_tree) +
  labs(fill = "SC") +
  xlab(element_blank()) + ylab("haplogroup proportions") +
  theme(text = element_text(size = 25)) +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "white")) +
  coord_flip()
SC_AFR_barplot

# Rishishwar frequency barplot
#change order of factors for matching across groupings
metainfo_afr$SCL_haplogroups <- factor(metainfo_afr$SCL_haplogroups, levels = c("A", "B", "C", "D", "H", "J", "M", "U", "L3", "L4", "L2", "L5", "L1", "L0")) 

#plot
SCL_AFR_barplot <- ggplot(metainfo_afr) + 
  geom_bar(mapping = aes(x = Superpopulation.code, fill = SCL_haplogroups, y = (..count..)/sum(..count..)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values=colors_SCL_tree) +
  labs(fill = "SC+L") +
  xlab(element_blank()) + ylab("haplogroup proportions") +
  theme(text = element_text(size = 25)) +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "white")
  ) +
  coord_flip()
SCL_AFR_barplot

# save Yamamoto + Rishishwar barplots together for figure layout
barplots_SC_SCL <- ggarrange(SC_AFR_barplot, SCL_AFR_barplot, 
                                 labels = c("A", "B"),
                                 ncol = 1, nrow = 2)
barplots_SC_SCL
ggsave("frequency_barplots_yama_rish_afr.pdf", barplots_SC_SCL, limitsize = FALSE)


# rhierBAPS level1 frequency barplot
barplots_rhb1 <- ggplot(metainfo_afr) + 
  geom_bar(mapping = aes(x = Superpopulation.code, fill = factor(rhierBAPS_lev1), y = (..count..)/sum(..count..)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values=colors_rhb1_tree) +
  labs(fill = "rhB level 1") +
  xlab(element_blank()) + ylab("haplogroup proportions") +
  theme(text = element_text(size = 25)) +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "white")
  ) +
  coord_flip()
barplots_rhb1

# save rhierBAPS level 1 in same way for figure layout
barplots_rhb1 <- ggarrange(barplots_rhb1, 
                                labels = c("C", " "),
                                ncol = 1, nrow = 2)
barplots_rhb1
ggsave("frequency_barplots_rhb_level1_afr.pdf", barplots_rhb1, limitsize = FALSE)

# rhierBAPS level2 frequency barplot
#change order of factors for matching across groupings
metainfo_afr$rhierBAPS_lev2 <- factor(metainfo_afr$rhierBAPS_lev2, levels = c("1", "2", "3", "6", "4", "5", "7")) 

#plot
barplots_rhb2 <- ggplot(metainfo_afr) + 
  geom_bar(mapping = aes(x = Superpopulation.code, fill = factor(rhierBAPS_lev2), y = (..count..)/sum(..count..)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values=colors_rhb2_tree) +
  labs(fill = "rhB level 2") +
  xlab(element_blank()) + ylab("haplogroup proportions") +
  theme(text = element_text(size = 25)) +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "white")
  ) +
  coord_flip()
barplots_rhb2

# rhierBAPS level3 frequency barplot
#change order of factors for matching across groupings
metainfo_afr$rhierBAPS_lev3 <- factor(metainfo_afr$rhierBAPS_lev3, levels = c("1", "2", "4", "3", "13", "12", "11", "10", "6", "5", "9", "8", "7", "14", "15", "16")) 

#plot
barplots_rhb3 <- ggplot(metainfo_afr) + 
  geom_bar(mapping = aes(x = Superpopulation.code, fill = factor(rhierBAPS_lev3), y = (..count..)/sum(..count..)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values=colors_rhb3_tree) +
  labs(fill = "rhB level 3") +
  xlab(element_blank()) + ylab("haplogroup proportions") +
  theme(text = element_text(size = 25)) +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "white")
  ) +
  coord_flip()
barplots_rhb3

# save rhb2+3 barplots together for figure layout
barplots_rhb2_3 <- ggarrange(barplots_rhb2, barplots_rhb3, 
                                labels = c("D", "E"),
                                ncol = 1, nrow = 2)
barplots_rhb2_3
ggsave("frequency_barplots_rhb_level2-3_afr.pdf", barplots_rhb2_3, limitsize = FALSE)
