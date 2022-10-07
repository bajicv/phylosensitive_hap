## this script:
# creates a metainfo file for the mtDNA from 1000 Genomes Project
# adds haplogroup information acquired by running HaploGrep2 on mtDNA 
# removes known related individuals

## this script requires: 
# the following R packages:
library(tidyverse) #data handling

# loading the following data from 1KGP
IGSR_samples <- read_tsv(file = "igsr_samples_all.tsv") #https://www.internationalgenome.org/data-portal/sample contains all data, not just mtDNA
haplogroups_mtDNA <- read.table("1KGP_mtDNA_nogaps_haplogrep.txt", dec = ".", sep = "\t", header =TRUE) #from running HaploGrep2 on 1KGP_mtDNA_nogaps.fasta 
relatives <- read.table("20140625_related_individuals.txt", dec = ".", sep = "\t", header =TRUE) #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

#################################################

# merge haplogroup information with metainfo by accession numbers
metainfo <- left_join(haplogroups_mtDNA, IGSR_samples, by = c("SampleID" = "Sample name")) #automatically excludes every sample from IGSR that is not in mtDNA dataset

# remove relatives
relatives_IDs <- relatives$Sample # make string of accession numbers to be removed
metainfo <- metainfo[!grepl(paste(relatives_IDs, collapse='|'), metainfo$SampleID),] # remove from metainfo 

# create metainfo file
write.table(metainfo, file = "metainfo_OG.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
