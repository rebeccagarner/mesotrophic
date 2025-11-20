# Assign taxonomy with custom Silva_curated+Ramos+BOLD Cyanobacteria database
# https://benjjneb.github.io/dada2/bigdata_paired.html

# Load packages
library(dada2)
packageVersion("dada2")

# Load sequence table
load("cyano_sediment_seqtab_nochim.rda")
dim(seqtab_nochim)

# Assign taxonomy
taxonomy <- assignTaxonomy(seqtab_nochim, "silva138_ramos_bold_rmduplicates.fasta.gz", multithread = TRUE, minBoot = 80, taxLevels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
save(taxonomy, file = "cyano_sediment_taxonomy_silva138curated_ramos_bold_minbt80.rda")

# Signal end of program
print("Taxonomic assignment complete!")
