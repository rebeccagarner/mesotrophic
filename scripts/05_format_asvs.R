# Format ASVs

setwd("~/Google Drive/paul_and_rebeccas_projects/2022 LP sediment project/")

# Load libraries
library(tidyverse)
library(seqinr)


#### Load and format ASV data ####
# Load sequence table
load("output/dada2/cyano_sediment_seqtab_nochim.rda")

# Load taxonomy
load("output/dada2/cyano_sediment_taxonomy_silva138curated_ramos_bold_minbt80.rda")

# Format taxonomy table
taxonomy <- taxonomy %>%
  as_tibble(rownames = "sequence")
taxonomy[taxonomy == "NA"] <- NA


#### Assign ASV codes ####
# Number of unique ASVs
(nasvs <- n_distinct(taxonomy$sequence))

# Assign ASV codes to unique sequences
taxonomy <- taxonomy %>%
  mutate(asv_code = paste0("ASV", str_pad(string = 1:nasvs, width = nchar(nasvs), side = "left", pad = "0")))

# Write ASVs to fasta file
# write.fasta(sequences = as.list(taxonomy$sequence),
#             names = taxonomy$asv_code,
#             file.out = "output/dada2/cyano_sediment_seqtab_nochim.fasta")


#### Combine sequence counts and taxonomy ####
# Convert sequence table to long format and join taxonomy
asv_melt <- seqtab_nochim %>%
  as_tibble(rownames = "sample_id") %>%
  pivot_longer(!sample_id, names_to = "sequence", values_to = "nseqs") %>%
  filter(nseqs > 0) %>%
  left_join(taxonomy, by = "sequence")
sum(asv_melt$nseqs) == sum(seqtab_nochim)  # Should evaluate to TRUE

# Parse sample IDs
asv_melt <- asv_melt %>%
  mutate(sample_id = str_remove(sample_id, "^lp2017_")) %>%
  mutate(sample_id = str_remove(sample_id, "_16S$")) %>%
  separate(sample_id, into = c("lake_id1", "lake_id2", "midpoint1", "midpoint2"), sep = "-") %>%
  unite("lake_id", lake_id1, lake_id2, sep = "-") %>%
  unite("midpoint", midpoint1, midpoint2, sep = ".") %>%
  mutate(sample_id = str_c(lake_id, midpoint, sep = "_")) %>%
  select(sample_id, lake_id, midpoint, asv_code, sequence, kingdom:genus, nseqs)

# Write melted sequence data to file
# asv_melt %>%
#   write_tsv("output/dada2/cyano_sediment_melt_all.tsv", col_names = TRUE)
