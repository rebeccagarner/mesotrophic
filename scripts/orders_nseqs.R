setwd("~/reb.garner@gmail.com - Google Drive/My Drive/projects/paul_and_rebeccas_projects/2022 LP sediment project/")

source("scripts/08_site_by_genus.R")

# Chroococcales, Nostocales, Synechococcales, and Oscillatoriales
orders_target <- c("Chroococcales", "Nostocales", "Synechococcales", "Oscillatoriales")

orders_long <- site_by_order %>%
  pivot_longer(!c(sample_id, lake_id, midpoint), names_to = "order", values_to = "nseqs")

orders_long %>%
  distinct(order) %>%
  filter(order %in% orders_target)

orders_long %>%
  group_by(order) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  arrange(-nseqs)

orders_long %>%
  group_by(order) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  arrange(-nseqs) %>%
  filter(order %in% orders_target) %>%
  summarize(nseqs = sum(nseqs))

