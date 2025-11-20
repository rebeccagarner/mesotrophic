# Taxonomic composition

setwd("~/Google Drive/paul_and_rebeccas_projects/2022 LP sediment project/")

# Load libraries
library(tidyverse)
library(vegan)

# Load palettes
#source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted ASV table
cyanobacteriia_melt <- read_tsv("output/cyanobacteriia/cyano_sediment_melt_cyanobacteriia.tsv", col_names = TRUE)

# Format taxonomy
taxonomy <- cyanobacteriia_melt %>%
  distinct(asv_code, sequence, kingdom, phylum, class, order, family, genus)


#### Visualize taxonomic composition ####
# Calculate relative sequence abundance
samples_nseqs <- cyanobacteriia_melt %>%
  group_by(sample_id) %>%
  summarize(sample_nseqs = sum(nseqs)) %>%
  ungroup()

cyanobacteriia_melt <- cyanobacteriia_melt %>%
  left_join(samples_nseqs, by = "sample_id") %>%
  mutate(relseqs = nseqs/sample_nseqs * 100)
sum(cyanobacteriia_melt$relseqs)/100 == length(unique(cyanobacteriia_melt$sample_id))  # Should evaluate to TRUE

# Define colour palette for Cyanobacteriia families


# Sediment stratigraphic profiles
(taxonomy_barplots <- cyanobacteriia_melt %>%
    mutate(taxon = case_when(!is.na(family) ~ family,
                             is.na(family) & !is.na(order) ~ str_c("Unassigned ", order),
                             is.na(order) & !is.na(class) ~ str_c("Unassigned Cyanobacteriia family"))) %>%
    group_by(lake_id, midpoint, taxon) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ggplot() +
    facet_wrap(~lake_id,
               #scales = "free",
               nrow = 1) +
    geom_bar(aes(x = -midpoint, y = relseqs,
                 fill = taxon),
             stat = "identity") +
    #scale_fill_manual(values = palette_division, na.value = "black") +
    coord_flip() +  # coord_flip(expand = 0)
    labs(x = "Sediment depth (cm)",
         y = "Relative sequence abundance",
         fill = "Taxonomy") +
    theme_bw())
#ggsave("figures/cyano_sediment_barplots.pdf", taxonomy_barplots, "pdf", width = 16, height = 10, units = "in")


#### PCA ####
asvtable <- cyanobacteriia_melt %>%
  select(sample_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("sample_id")

asvtable_hellinger <- decostand(asvtable, "hellinger")

pca <- rda(asvtable_hellinger)

pca_sites <- scores(pca)$sites %>%
  as_tibble(rownames = "sample_id") %>%
  separate(sample_id, into = c("lake_id", "midpoint"), sep = "_")
pca_species <- scores(pca)$species %>%
  as_tibble(rownames = "asv_code")
(pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1))
(pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1))

# PC contributions
contributors <- pca_species %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(magnitude, n = 10) %>%
  pull(asv_code)

# Curate taxon contributor labels in PCA plot
contribution <- pca_species %>%
  filter(asv_code %in% contributors) %>%
  left_join(taxonomy, by = "asv_code") %>%
  mutate(label = case_when(!is.na(genus) ~ str_c(asv_code, " ", genus, "\n(", order, ")", sep = ""),
                           is.na(genus) & !is.na(family) ~ str_c(asv_code, " ", family, "\n(", order, ")", sep = ""),
                           is.na(family) & !is.na(order) ~ str_c(asv_code, " ", order, sep = ""),
                           is.na(order) & !is.na(class) ~ str_c(asv_code, " ", class, sep = ""))) %>%
  mutate(label = str_replace_all(label, "_", " "))

# Plot PCA
(pca_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#E6E6E6") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#E6E6E6") +
    geom_text(aes(x = 1.1*PC1, y = 1.1*PC2, label = label),
              data = contribution,
              colour = "#A0A0A0", size = 4, fontface = "italic") +
    geom_segment(aes(x = 0, y = 0, xend = 1*PC1, yend = 1*PC2),
                 data = contribution,
                 arrow = arrow(length = unit(x = 0.2, units = "cm")),
                 alpha = 0.7, colour = "#A0A0A0") +
    labs(x = paste("PC axis 1 (", pc1_percent_var, "%)", sep = ""),
         y = paste("PC axis 2 (", pc2_percent_var, "%)", sep = ""),
         colour = "Lake") +
    geom_point(aes(x = PC1, y = PC2,
                   colour = lake_id),
               data = pca_sites, size = 5, alpha = 0.7, stroke = 2) +
    geom_text(aes(x = PC1, y = PC2, label = midpoint),
              data = pca_sites,
              colour = "black", size = 2, fontface = "bold", check_overlap = TRUE) +
    #scale_colour_manual(values = palette_lake) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          legend.key = element_blank(), legend.position = "right",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5)))#ggsave("figures/ela18s_pca.pdf", pca_plot, "pdf", width = 15, height = 12, units = "in")
#ggsave("figures/cyano_sediment_pca.pdf", pca_plot, "pdf", width = 12, height = 10, units = "in")
