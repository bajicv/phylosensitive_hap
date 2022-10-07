## this script:
# runs rhierBAPS with the settings described in the thesis
# saves the output for usage in downstream analyses

## this script requires: 
# the following R packages:
library(rhierbaps)

# loading FASTA (created by the appropriate creating_FASTA_subset_xyz.R) as a SNP matrix
snp_matrix <- load_fasta("KGP_fasta_african_subset.fasta", keep.singletons = T) 

#################################################

# run rhierBAPS 
rhb_results <- hierBAPS(snp_matrix, max.depth = 10, quiet = FALSE) # n.pop = default = number of individuals/5 (does not need specification)

# save partition dataframe from results 
write.csv(rhb_results$partition.df, file = "rhb_african_level10_singletons.csv", 
          col.names = TRUE, row.names = FALSE)