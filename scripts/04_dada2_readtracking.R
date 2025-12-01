# Read tracking through DADA2 pipeline

library(tidyverse)

# Import read tracking file
readtracking <- read_tsv("output/dada2/cyano_sediment_dada2_readtracking.tsv", col_names = TRUE)

load("output/dada2/cyano_sediment_seqtab_nochim.rda")
dim(seqtab_nochim)

addChim <- function(readtracking, seqtab_nochim) {
  readtracking %>%
    left_join(tibble(sample_id = rownames(seqtab_nochim),
                     nonchim = rowSums(seqtab_nochim)), by = "sample_id") %>%
    arrange(sample_id) %>%
    mutate(pct_seqs_retained = nonchim/filtered * 100) %>%
    return()
}

(readtracking <- addChim(readtracking, seqtab_nochim))

load("output/dada2/cyano_sediment_seqtab_merged.rda")
dim(seqtab)

calculateASVs <- function(seqtab_merged, seqtab_nochim) {
  seqtab_merged %>%
    as_tibble(rownames = "sample_id") %>%
    arrange(sample_id) %>%
    pivot_longer(!sample_id, names_to = "sequence", values_to = "nseqs") %>%
    filter(nseqs > 0) %>%
    mutate(seq_type = case_when(sequence %in% colnames(seqtab_nochim) ~ "nonchim",
                                TRUE ~ "chim")) %>%
    group_by(sample_id, seq_type) %>%
    count(name = "nasvs") %>%
    ungroup() %>%
    pivot_wider(names_from = seq_type, values_from = nasvs, names_glue = "{.value}_{seq_type}") %>%
    arrange(sample_id) %>%
    mutate(nasvs_total = nasvs_chim + nasvs_nonchim) %>%
    return()
}

(nasvs <- calculateASVs(seqtab, seqtab_nochim))

readtracking <- readtracking %>%
  left_join(nasvs, by = "sample_id")

# readtracking %>%
#   write_tsv("output/dada2/cyano_sediment_dada2_readtracking_all.tsv")
