
setwd("~/reb.garner@gmail.com - Google Drive/My Drive/projects/paul_and_rebeccas_projects/2022 LP sediment project/")

# Load libraries
library(tidyverse)


orders_target <- c("Chroococcales", "Nostocales", "Synechococcales", "Oscillatoriales")

asv_melt %>%
  filter(class == "Cyanobacteriia") %>%
  filter(order %in% orders_target) %>%
  group_by(order) %>%
  summarize(nasvs = length(unique(asv_code)),
            nseqs = sum(nseqs)) %>%
  ungroup() %>%
  arrange(-nseqs)

asv_melt %>%
  filter(class == "Cyanobacteriia") %>%
  group_by(sample_id) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  summarize(min_nseqs = min(nseqs),
            max_nseqs = max(nseqs),
            mean_nseqs = mean(nseqs),
            median_nseq = median(nseqs))


asv_melt %>%
  filter(class == "Cyanobacteriia") %>%
  mutate(order_class = case_when(order %in% orders_target ~ "orders4",
                                 TRUE ~ "other")) %>%
  group_by(sample_id, order_class) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  pivot_wider(names_from = order_class, values_from = nseqs, values_fill = 0)%>%
  mutate(total_seqs = orders4 + other) %>%
  mutate(pct_order4 = orders4/total_seqs*100) %>%
  summarize(min_pseqs = min(pct_order4),
            max_pseqs = max(pct_order4),
            mean_pseqs = mean(pct_order4),
            median_pseq = median(pct_order4))



