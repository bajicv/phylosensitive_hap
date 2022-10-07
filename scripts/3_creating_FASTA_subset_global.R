## this script:
# subsets a FASTA file based on random sampling from every population in the 1000 Genomes Project (1KGP) + RSRS for iqtree root
# creates a FASTA file used as input for various downstream analyses

## this script requires: 
# the following R packages:
library(tidyverse) #data handling
library(seqinr) #loading FASTA file into R

# loading the following data 1KGP
metainfo <- read.table("metainfo_OG.txt", header=T, sep=" ") #metainfo file (created in xyz.R)
KGP_fasta <- read.fasta("1KG_rCRS_RSRS_norelatives_nopolyc_NMreplaced_mafft_indelcorr_postmafft_new.fasta") #edited/filtered mtDNA (created in xyz.R), multiple sequence alignment with post-alignment indel correction (described in thesis methods)

#################################################

# check populations present for overview
unique(metainfo$Population.code) #contains 1 individual with mixed population designation "IBS,MSL"
# remove individual with mixed population designation, as no random sampling can be performed on this population
metainfo <-  filter(metainfo, metainfo$Population.code != "IBS,MSL") #remove mixed pop individual

# randomly sample 20 individuals / population from metainfo 
# reason: needed to subset the FASTA, as the FASTA only contains sample IDs and no population info
all_subset <- data.frame(matrix(ncol = 0, nrow = 0))
for (i in unique(metainfo$Population.code)) {
  all_subset <- rbind(all_subset, sample_n((subset(metainfo, Population.code == i)), 20))
}

# subset FASTA based on sample IDs from random sampling
KGP_fasta_all_subset <- KGP_fasta[names(KGP_fasta) %in% all_subset$SampleID] 

# create FASTA as input for rhierBAPS
write.fasta(sequences = KGP_fasta_all_subset, names = names(KGP_fasta_all_subset), 
            file.out = "KGP_fasta_all_subset.fasta") #input for rhierBAPS

# extract RSRS from original FASTA and add it to the FASTA containing randomly sampled individuals
# reason: RSRS is used to root the tree in iqtree, but is not needed for rhierBAPS
KGP_fasta_RSRS <- KGP_fasta[names(KGP_fasta) == "RSRS"] 
KGP_fasta_RSRS_all <- c(KGP_fasta_RSRS, KGP_fasta_all_subset)

# create FASTA including RSRS as input for iqtree
write.fasta(sequences = KGP_fasta_RSRS_all, names = names(KGP_fasta_RSRS_all), 
            file.out = "KGP_fasta_RSRS_all_iqtree.fasta") #input for iqtree 
