# Format ASVs

#setwd("~/Google Drive/paul_and_rebeccas_projects/2022 LP sediment project/")
setwd("~/reb.garner@gmail.com - Google Drive/My Drive/projects/paul_and_rebeccas_projects/2022 LP sediment project/")

# Load libraries
library(tidyverse)


#### Load and format data ####
# Import melted sequence data to file
asv_melt <- read_tsv("output/dada2/cyano_sediment_melt_all.tsv", col_names = TRUE)
length(unique(asv_melt$sample_id))  # Number of samples (86)
length(unique(asv_melt$asv_code))  # Number of ASVs (14588)
sum(asv_melt$nseqs)  # Number of total sequences (4942025)


#### Visualize data set quality ####
# Plot histogram of n sequences assigned to each ASV
asv_melt %>%
  group_by(asv_code, sequence) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  ggplot() +
  geom_histogram(aes(x = nseqs, y = ..count..), binwidth = 100)

# Plot pie chart of n sequences per ASV
(asvs_nseqs_piechart <- asv_melt %>%
    group_by(asv_code, sequence) %>%
    summarize(nseqs = sum(nseqs)) %>%
    ungroup() %>%
    mutate(sequence_nseqs = case_when(nseqs == 1 ~ "singleton",
                                      nseqs == 2 ~ "doubleton",
                                      nseqs > 2 & nseqs <= 10 ~ "3 - 10 seqs",
                                      nseqs > 10 & nseqs <= 100 ~ "11 - 100 seqs",
                                      nseqs >= 100 ~ ">100 seqs")) %>%
    group_by(sequence_nseqs) %>%
    tally(name = "nasvs") %>%
    ggplot(aes(x = "", y = nasvs, fill = sequence_nseqs)) +
    geom_bar(stat = "identity") +
    coord_polar("y") +
    geom_text(aes(label = str_c(nasvs, "\n(", round(nasvs/(length(unique(asv_melt$asv_code)))*100, 0), "%)"), x = 1),
              position = position_stack(vjust = 0.5),
              size = 4) +
    theme_void())
#ggsave("figures/cyano_sediment_asv_nseqs_piechart.pdf", asvs_nseqs_piechart, "pdf", width = 5, height = 4, units = "in")

# Plot histogram of sequence counts by sample
(samples_nseqs_histogram <- asv_melt %>%
    group_by(sample_id) %>%
    summarize(nseqs = sum(nseqs)) %>%
    ungroup() %>%
    ggplot() +
    geom_histogram(aes(x = nseqs, y = ..count..), binwidth = 1000) +
    labs(x = "Number of sequences",
         y = "Number of samples") +
    theme_bw())
#ggsave("figures/cyano_sediment_samples_nseqs_histogram.pdf", samples_nseqs_histogram, "pdf", width = 5, height = 4, units = "in")


#### Curate ASVs ####
# Remove global singleton and doubleton ASVs
asvs_rm <- asv_melt %>%
  group_by(asv_code, sequence) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  filter(nseqs <= 2) %>%
  pull(sequence)

asv_melt %>%
  filter(sequence %in% asvs_rm) %>%
  summarize(nseqs = sum(nseqs))  # Number of sequences assigned to singletons or doubletons (1962)

asv_melt <- asv_melt %>%
  filter(!sequence %in% asvs_rm)
length(unique(asv_melt$sample_id))  # Number of samples (86)
length(unique(asv_melt$asv_code))  # Number of retained ASVs (13508)
sum(asv_melt$nseqs)  # Number of total sequences (4940063)


# #### Curate samples ####
# # Remove samples with fewer than 10,000 sequences
# samples_rm <- asv_melt %>%
#   group_by(sample_id) %>%
#   summarize(nseqs = sum(nseqs)) %>%
#   ungroup() %>%
#   filter(nseqs < 10000) %>%
#   pull(sample_id)
# samples_rm  # Remove samples "06-199_8.25" and "06-199_9.25"
# 
# asv_melt <- asv_melt %>%
#   filter(!sample_id %in% samples_rm)
# length(unique(asv_melt$sample_id))  # Number of retained samples


#### Curate ASVs based on taxonomy ####
# Create taxonomy table
taxonomy <- asv_melt %>%
  select(asv_code, sequence, kingdom:genus) %>%
  distinct(asv_code, .keep_all = TRUE)

taxonomy %>%
  pivot_longer(!c(asv_code, sequence), names_to = "rank", values_to = "taxon") %>%
  mutate(assignment = case_when(is.na(taxon) ~ "unassigned",
                                !is.na(taxon) ~ "assigned")) %>%
  group_by(rank, assignment) %>%
  count(name = "nasvs") %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = factor(rank, levels = c("kingdom", "phylum", "class", "order", "family", "genus")),
               y = nasvs,
               fill = factor(assignment, levels = c("unassigned", "assigned"))),
           stat = "identity") +
  labs(fill = "Assignment") +
  theme(axis.title.x = element_blank())

# Plot pie chart of all taxa in dataset
(phyla_nseqs_piechart <- asv_melt %>%
  mutate(taxgroup = case_when(order == "Chloroplast" ~ "Chloroplast",
                              class == "Sericytochromatia" | class == "Vampirivibrionia" ~ "Non-photosynthetic Cyanobacteria",
                              class == "Cyanobacteriia" ~ "Cyanobacteriia",
                              phylum == "Cyanobacteria" ~ "Unassigned Cyanobacteria",
                              is.na(phylum) ~ "Unassigned phylum",
                              TRUE ~ "Other phylum")) %>%
  group_by(taxgroup) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = nseqs,
             fill = taxgroup)) +
  geom_bar(stat = "identity") +
  coord_polar("y") +
  geom_text(aes(label = str_c(nseqs, "\n(", round(nseqs/sum(asv_melt$nseqs)*100, 0), "%)"), x = 1.2),
            position = position_stack(vjust = 0.5),
            size = 3) +
  labs(fill = "Taxonomy") +
  theme_void())
#ggsave("figures/cyano_sediment_phyla_nseqs_piechart.pdf", phyla_nseqs_piechart, "pdf", width = 5, height = 4, units = "in")

# Identify ASVs assigned to Chloroplasts
chloroplasts <- taxonomy %>%
  filter(order == "Chloroplast") %>%
  pull(sequence)
length(unique(chloroplasts))  # Number of ASVs assigned to Chloroplasts (383)

asv_melt %>%
  filter(sequence %in% chloroplasts) %>%
  summarize(nseqs = sum(nseqs))  # Number of sequences assigned to Chloroplasts (582874)

# Examine ASVs assigned to non-photosynthetic Cyanobacteria
taxonomy %>%
  filter(class == "Sericytochromatia" | class == "Vampirivibrionia") %>%
  nrow()  # Number of ASVs assigned to non-photosynthetic Cyanobacteria (130)

asv_melt %>%
  filter(class == "Sericytochromatia" | class == "Vampirivibrionia") %>%
  summarize(nseqs = sum(nseqs))  # Number of sequences assigned to non-photosynthetic Cyanobacteria (11644)

# Filter Cyanobacteriia  (remove chloroplasts)
cyanobacteriia_melt <- asv_melt %>%
  filter(class == "Cyanobacteriia") %>%
  filter(!sequence %in% chloroplasts)
length(unique(cyanobacteriia_melt$sample_id))  # Number of samples (82)
length(unique(cyanobacteriia_melt$asv_code))  # Number of ASVs (882)
sum(cyanobacteriia_melt$nseqs)  # Number of sequences (789931)

# Write melted sequence table to file
# cyanobacteriia_melt %>%
#   write_tsv("output/cyanobacteriia/cyano_sediment_melt_cyanobacteriia.tsv", col_names = TRUE)

# Write sequence table (site by ASV) to file
# cyanobacteriia_melt %>%
#   select(sample_id, lake_id, midpoint, asv_code, nseqs) %>%
#   pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
#   write_tsv("output/cyanobacteriia/cyano_sediment_asvtable_cyanobacteriia.tsv", col_names = TRUE)

# Write taxonomy table to file
# cyanobacteriia_melt %>%
#   select(asv_code, sequence, kingdom:genus) %>%
#   distinct(asv_code, .keep_all = TRUE) %>%
#   write_tsv("output/cyanobacteriia/cyano_sediment_taxonomy_cyanobacteriia.tsv", col_names = TRUE)

# Which samples have been retained?
samples_retained <- cyanobacteriia_melt %>%
  distinct(sample_id, lake_id, midpoint)

# Which samples have been thrown out as a result of taxonomy curation?
asv_melt %>%
  distinct(sample_id) %>%
  filter(!sample_id %in% samples_retained$sample_id)
# 06-199_8.25, 08-160_22.5, 08-160_25.5, 17-060_23.5
