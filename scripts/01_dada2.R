# Run DADA2 pipeline on Cutadapt-trimmed reads
# For LakePulse 2019 surface water 16S rRNA gene amplicon sequencing data
# https://benjjneb.github.io/dada2/bigdata_paired.html

setwd("~/Google Drive/paul_and_rebeccas_projects/2022 LP sediment project/filtered_samples/")

# Load packages
library(dada2)
packageVersion("dada2")
library(tidyverse)

# Print system time
print(paste0("Program start: ", Sys.time()))

# Identify 2019 samples
(samples <- scan("sediment2022_samples.txt", what = "character"))

# Locate and name DADA2-filtered forward and reverse read files
(filtFs <- tibble(file_name = sort(list.files(pattern = "_sub_R1_filtered.fq.gz"))) %>%
    filter(file_name %in% paste0(samples, "_sub_R1_filtered.fq.gz")) %>%
    pull(file_name))
names(filtFs) <- samples

(filtRs <- tibble(file_name = sort(list.files(pattern = "_sub_R2_filtered.fq.gz"))) %>%
    filter(file_name %in% paste0(samples, "_sub_R2_filtered.fq.gz")) %>%
    pull(file_name))
names(filtRs) <- samples

# Learn error rates
print("Learning error rates...")
set.seed(33)
err_forward_reads <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE)
err_reverse_reads <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# Plot error rates
pdf("cyano_sediment_errors_fwd.pdf", width = 10, height = 10)
plotErrors(err_forward_reads, nominalQ = TRUE)
dev.off()

pdf("cyano_sediment_errors_rev.pdf", width = 10, height = 10)
plotErrors(err_reverse_reads, nominalQ = TRUE)
dev.off()

# Track read loss through the pipeline
getN <- function(x) sum(getUniques(x))

summary_table <- data.frame(row.names = samples,
                            filtered = sapply(filtFs, getN),
                            dada_fwd = NA,
                            dada_rev = NA,
                            merged = NA)

# Dereplicate reads, infer ASVs (no sample pooling), and merge paired-end reads
print("Dereplicating reads, inferring ASVs, and merging paired-end reads...")
mergers <- vector("list", length(samples))
names(mergers) <- samples
for(i in samples) {
  cat("Processing:", i, "\n")
  derepF <- derepFastq(filtFs[[i]], verbose = TRUE)
  ddF <- dada(derepF, err = err_forward_reads, pool = FALSE, multithread = TRUE)
  summary_table$dada_fwd[which(rownames(summary_table) == i)] <- getN(ddF)
  
  derepR <- derepFastq(filtRs[[i]], verbose = TRUE)
  ddR <- dada(derepR, err = err_reverse_reads, pool = FALSE, multithread = TRUE)
  summary_table$dada_rev[which(rownames(summary_table) == i)] <- getN(ddR)
  
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE, minOverlap = 20)
  mergers[[i]] <- merger
  summary_table$merged[which(rownames(summary_table) == i)] <- getN(merger)
}

# Construct sequence table
print("Constructing sequence table...")
seqtab <- makeSequenceTable(mergers)
save(seqtab, file = "cyano_sediment_seqtab_merged.rda")
dim(seqtab)

# Save read tracking table
summary_table %>%
  as_tibble(rownames = "sample_id") %>%
  write_tsv("cyano_sediment_dada2_readtracking.tsv", col_names = TRUE)

# Print system time
print(paste0("Program end: ", Sys.time()))
