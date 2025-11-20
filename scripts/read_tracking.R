setwd("~/reb.garner@gmail.com - Google Drive/My Drive/projects/paul_and_rebeccas_projects/2022 LP sediment project/")

library(tidyverse)
library(readxl)
library(janitor)

readtracking <- read_tsv("output/dada2/cyano_sediment_dada2_readtracking_all.tsv")
genomequebec <- read_xlsx("data/LakePulse_-_Cyano16S_COI_amplicon_sediments_MiSeqReadSet_2022-03-04.xlsx") %>%
  clean_names() %>%
  mutate(sample_id = str_replace_all(name, "\\.", "-")) %>%
  filter(grepl("16S", sample_id)) %>%
  mutate(raw = as.numeric(str_remove_all(number_of_reads, ",")))

data <- readtracking %>%
  left_join(genomequebec, by = c("sample_id"))

data %>%
  mutate(filt_raw_pct = filtered/raw*100) %>%
  select(filt_raw_pct) %>%
  summarize(min = min(filt_raw_pct),
            max = max(filt_raw_pct),
            mean = mean(filt_raw_pct),
            median = median(filt_raw_pct))

