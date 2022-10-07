## this script:
# subsets a FASTA file based on the superpopulation "AFR" from the 1000 Genomes Project (1KGP) + RSRS for iqtree root
# creates a FASTA file used as input for various downstream analyses

## this script requires: 
# the following R packages:
library(tidyverse) #data handling
library(seqinr) #loading FASTA file into R

# loading the following data from 1KGP
metainfo <- read.table("metainfo_OG.txt", header=T, sep=" ") #metainfo file (created in metainfo.R)
KGP_fasta_raw <- read.fasta("1KG_rCRS_RSRS_norelatives_nopolyc_NMreplaced_mafft_indelcorr_postmafft_new.fasta") #edited/filtered mtDNA (created in xyz.R), multiple sequence alignment with post-alignment indel correction (described in thesis methods)

#################################################

# check populations for overview
unique(metainfo$Superpopulation.code)

# subset metainfo by superpopulation == AFR to acquire sample IDs assigned AFR
# reason: needed to subset the FASTA, as the FASTA only contains sample IDs and no superpopulation info
metainfo_afr <- metainfo %>% 
  subset(Superpopulation.code == "AFR") 

# subset FASTA based on sample IDs from african metainfo
KGP_fasta_afr_subset <- KGP_fasta_raw[names(KGP_fasta_raw) %in% metainfo_afr$SampleID] 

# create FASTA as input for rhierBAPS
write.fasta(sequences = KGP_fasta_afr_subset, names = names(KGP_fasta_afr_subset), 
            file.out = "KGP_fasta_african_subset.fasta") 

# extract RSRS from original FASTA and add it to the AFR subset FASTA
# reason: RSRS is used to root the tree in iqtree, but is not needed for rhierBAPS
KGP_fasta_RSRS <- KGP_fasta_raw[names(KGP_fasta_raw) == "RSRS"] 
KGP_fasta_RSRS_afr <- c(KGP_fasta_RSRS, KGP_fasta_afr_subset) 

# create FASTA including RSRS as input for iqtree
write.fasta(sequences = KGP_fasta_RSRS_afr, names = names(KGP_fasta_RSRS_afr), 
            file.out = "KGP_fasta_RSRS_african_iqtree.fasta") #input for iqtree 




