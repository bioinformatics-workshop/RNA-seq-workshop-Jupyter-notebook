#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
library(glue)

# Read in file
data <- read_tsv("raw/PRJNA950346.metadata.tmp", col_names = F)

# Clean up file
data_clean <- data %>%
       select(srr_id = X8, 
              ecotype = X3, 
              genotype = X4,
              treatment = X5,
              tissue = X2) %>%
       mutate(biorep = c(rep(1:3,2))) %>%
       mutate(genotype = case_when(
              grepl("WT", genotype) ~ "wt",
              TRUE ~ "miR163_mut"
              )) %>%
       mutate(samplename = case_when(
              grepl("wt", genotype) ~ glue("wt_{biorep}"),
              TRUE ~ glue("mir163_{biorep}")
              )) %>%
       mutate(fq1 = glue("raw/{srr_id}_1.fastq.gz"),
              fq2 = glue("raw/{srr_id}_2.fastq.gz")) %>%
       select(samplename, fq1, fq2, srr_id, ecotype, genotype, treatment, tissue, biorep)

# Export file
write_csv(data_clean, "metadata/metadata.csv")

### Metadata for the subsampled dataset

# Clean up file
data_clean_sub1M <- data %>%
       select(srr_id = X8, 
              ecotype = X3, 
              genotype = X4,
              treatment = X5,
              tissue = X2) %>%
       mutate(srr_id = glue("{srr_id}_sub1M")) %>%
       mutate(biorep = c(rep(1:3,2))) %>%
       mutate(genotype = case_when(
              grepl("WT", genotype) ~ "wt", 
              TRUE ~ "miR163_mut"
              )) %>%
       mutate(samplename = case_when(
              grepl("wt", genotype) ~ glue("wt_{biorep}"),
              TRUE ~ glue("mir163_{biorep}")
              )) %>%
       mutate(fq1 = glue("raw/{srr_id}_1.fastq.gz"),
              fq2 = glue("raw/{srr_id}_2.fastq.gz")) %>%
       select(samplename, fq1, fq2, srr_id, ecotype, genotype, treatment, tissue, biorep)

# Export file
write_csv(data_clean_sub1M, "metadata/metadata_sub1M.csv")
