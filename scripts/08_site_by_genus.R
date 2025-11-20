# Format ASVs

#setwd("~/Google Drive/paul_and_rebeccas_projects/2022 LP sediment project/")
setwd("~/reb.garner@gmail.com - Google Drive/My Drive/projects/paul_and_rebeccas_projects/2022 LP sediment project/")

# Load libraries
library(tidyverse)


#### Load and format data ####
# Import melted ASV data
asv_melt <- read_tsv("output/cyanobacteriia/cyano_sediment_melt_cyanobacteriia.tsv", col_names = TRUE)
length(unique(asv_melt$sample_id))  # Number of samples
length(unique(asv_melt$asv_code))  # Number of ASVs
sum(asv_melt$nseqs)  # Number of total sequences


#### Format site by genus table(s) ####
# Extract genus when there is extra species information in the genus assignment
asv_melt <- asv_melt %>%
  rename(genus_species = genus) %>%
  mutate(genus_species = case_when(genus_species == "SU2_symbiont_group" ~ "SU2-symbiont-group",
                                   TRUE ~ genus_species)) %>%
  mutate(genus = str_remove_all(genus_species, "_.*"))

# Collapse sequences by sample and genus
genus_melt <- asv_melt %>%
  group_by(sample_id, lake_id, midpoint, kingdom, phylum, class, order, family, genus, genus_species) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup()

# Rename NA genus to "Unclassified_genus"
genus_melt <- genus_melt %>%
  mutate(genus = case_when(is.na(genus) ~ "Unclassified_genus",
                           TRUE ~ genus),
         genus_species = case_when(is.na(genus_species) ~ "Unclassified_genus_species",
                                   TRUE ~ genus_species))

# Inspect ASVs with unclassified genera
genus_melt %>%
  group_by(genus) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  arrange(-nseqs) %>%
  mutate(pct_seqs = nseqs/sum(asv_melt$nseqs) * 100)

# Create site by genus_species matrix (all NA genus_species to a single "Unclassified_genus_species" column)
site_by_genus_species <- genus_melt %>%
  select(sample_id, lake_id, midpoint, genus_species, nseqs) %>%
  group_by(sample_id, lake_id, midpoint, genus_species) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  pivot_wider(names_from = genus_species, values_from = nseqs, values_fill = 0)
sum(site_by_genus_species[,-c(1:3)]) == sum(asv_melt$nseqs)  # Should evaluate to TRUE

# Write site by genus_species matrix to file
# site_by_genus_species %>%
#     write_tsv("output/cyanobacteriia/cyano_sediment_site_by_genus_species_cyanobacteriia.tsv", col_names = TRUE)

# Create site by genus matrix (all NA genera to a single "Unclassified_genus" column)
site_by_genus <- genus_melt %>%
  select(sample_id, lake_id, midpoint, genus, nseqs) %>%
  group_by(sample_id, lake_id, midpoint, genus) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  pivot_wider(names_from = genus, values_from = nseqs, values_fill = 0)
sum(site_by_genus[,-c(1:3)]) == sum(asv_melt$nseqs)  # Should evaluate to TRUE

# Write site by genus matrix to file
# site_by_genus %>%
#     write_tsv("output/cyanobacteriia/cyano_sediment_site_by_genus_cyanobacteriia.tsv", col_names = TRUE)

# Create site by genus matrix (NA genera to separate order-rank "Unclassified" columns for
# Chroococcales, Nostocales, Synechococcales, or Oscillatoriales)
site_by_genus_separate_unclassifieds <- genus_melt %>%
  mutate(genus = case_when(genus == "Unclassified_genus" & order == "Chroococcales" ~ "Unclassified_Chroococcales",
                           genus == "Unclassified_genus" & order == "Nostocales" ~ "Unclassified_Nostocales",
                           genus == "Unclassified_genus" & order == "Synechococcales" ~ "Unclassified_Synechococcales",
                           genus == "Unclassified_genus" & order == "Oscillatoriales" ~ "Unclassified_Oscillatoriales",
                           TRUE ~ genus)) %>%
  select(sample_id, lake_id, midpoint, genus, nseqs) %>%
  group_by(sample_id, lake_id, midpoint, genus) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  pivot_wider(names_from = genus, values_from = nseqs, values_fill = 0)
sum(site_by_genus_separate_unclassifieds[,-c(1:3)]) == sum(asv_melt$nseqs)  # Should evaluate to TRUE

# Write site by genus matrix with separate order-rank unclassified genera to file
# site_by_genus_separate_unclassifieds %>%
#     write_tsv("output/cyanobacteriia/cyano_sediment_site_by_genus_separate_unclassifieds_cyanobacteriia.tsv", col_names = TRUE)


#### Format site by order table(s) ####

# Collapse sequences by sample and order
order_melt <- asv_melt %>%
  group_by(sample_id, lake_id, midpoint, kingdom, phylum, class, order) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup()

# Rename NA order to "Unclassified_order"
order_melt <- order_melt %>%
  mutate(order = case_when(is.na(order) ~ "Unclassified_order",
                           TRUE ~ order))

# Inspect ASVs with unclassified orders
order_melt %>%
  group_by(order) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  arrange(-nseqs) %>%
  mutate(pct_seqs = nseqs/sum(asv_melt$nseqs) * 100)

# Create site by order matrix (all NA order to a single "Unclassified_order" column)
site_by_order <- order_melt %>%
  select(sample_id, lake_id, midpoint, order, nseqs) %>%
  group_by(sample_id, lake_id, midpoint, order) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  pivot_wider(names_from = order, values_from = nseqs, values_fill = 0)
sum(site_by_order[,-c(1:3)]) == sum(asv_melt$nseqs)  # Should evaluate to TRUE

# Write site by order matrix to file
# site_by_order %>%
#     write_tsv("output/cyanobacteriia/cyano_sediment_site_by_order_cyanobacteriia.tsv", col_names = TRUE)
