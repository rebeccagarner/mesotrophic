# Remove bimeras
# https://benjjneb.github.io/dada2/bigdata_paired.html

# Load packages
library(dada2)
packageVersion("dada2")

# Remove chimeras
load("cyano_sediment_seqtab_merged.rda")
dim(seqtab)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "pooled", multithread = TRUE)
save(seqtab_nochim, file = "cyano_sediment_seqtab_nochim.rda")
dim(seqtab_nochim)
