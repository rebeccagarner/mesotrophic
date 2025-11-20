setwd("~/Library/CloudStorage/Dropbox-Personal/projects/lpmags/")

# Load libraries
library(tidyverse)

metadata <- read_csv("data/environmental/lakepulse_data_curated_2021-05-04.csv", col_names = TRUE)

# Green Lake
# Lac de Saint-Damase
# Lac du Marin-à-Gouin
# Lac Michaud
# Lillabelle Lake
# Stoco Lake
lakes <- c("06-199",
           "08-160",
           "08-192",
           "17-059",
           "17-060",
           "17-075")

metadata <- metadata %>%
  filter(lakepulse_id %in% lakes)

