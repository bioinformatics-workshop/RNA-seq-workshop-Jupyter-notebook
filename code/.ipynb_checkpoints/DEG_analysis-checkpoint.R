#!/usr/bin/env Rscript

# Load helper functions
source("code/utility_scripts.R")

# load libraries
library(tidyverse)
library(DESeq2)
library(glue)
library(UpSetR)

# Load annotation file
annot_files <- c("genome/MEAM1_annotation_v1.2.xlsx",
                "genome/Rickettsia_annotation_v1.xlsx")

annot <- lapply(annot_files, readxl::read_xlsx) %>% bind_rows

# Load sample info and metadata
sampleinfo <- read_csv("metadata/metadata.csv", col_names = TRUE)

# Create factors for various comparisons
# factor-genotype (M1, M2)
# factor-sex (MA, FE)
# factor-infection (IN, UN)
# factor-substrate (CP, SU)
# factor-infection_substrate
# factor-infection_genotype
# factor-infection_sex
# factor-genotype_sex
# factor-infection_genotype_sex
# factor-infection_genotype_substrate
# factor-infection_sex_substrate
# factor-genotype_sex_substrate

sampleinfo <- sampleinfo %>%
    dplyr::rename(infection = infection_status,
                substrate = feeding_substrate) %>%
    dplyr::rename_with(~str_c("fct_",.), infection:substrate) %>%
    dplyr::mutate(
        fct_infection = ifelse(grepl("un", fct_infection),"UN","IN"),
        fct_genotype = ifelse(grepl("1", fct_genotype), "M1", "M2"),
        fct_sex = ifelse(grepl("female", fct_sex), "FE", "MA"),
        fct_substrate = ifelse(grepl("cowpea", fct_substrate), "CP", "SU")    
        ) %>%
    dplyr::mutate(
        fct_infection_substrate = glue("{fct_infection}_{fct_substrate}"),
        fct_infection_genotype = glue("{fct_infection}_{fct_genotype}"),
        fct_infection_sex = glue("{fct_infection}_{fct_sex}"),
        fct_genotype_sex = glue("{fct_genotype}_{fct_sex}"),
        fct_infection_genotype_sex = glue("{fct_infection}_{fct_genotype}_{fct_sex}"),
        fct_infection_genotype_substrate = glue("{fct_infection}_{fct_genotype}_{fct_substrate}"),
        fct_infection_sex_substrate = glue("{fct_infection}_{fct_sex}_{fct_substrate}"),
        fct_genotype_sex_substrate = glue("{fct_genotype}_{fct_sex}_{fct_substrate}"),
        fct_all = glue("{fct_infection}_{fct_genotype}_{fct_sex}_{fct_substrate}"),
        ) %>%
    column_to_rownames(var = "sample")

sample_order <- rownames(sampleinfo)

# Load count matrix file from featurecounts and clean up the header
fcnts <- read_tsv("analysis/featurecounts/allsamples_merged.fcnts.txt")

counts <- fcnts %>%
    column_to_rownames(var = "Geneid")

counts <- counts[,sample_order]

# Obtain deseq2 output per factor comparisons
deseq_per_factor <- function(fct){

    print(glue("Processing factor: {fct}"))

    print("generating comparisons")
    # get a list of all possible combinations from sample groupings
    cmp <- gtools::combinations(n = length(unique(sampleinfo[,fct])),
            r = 2,
            v = sampleinfo[,fct])
    cmp <- as.list(data.frame(t(cmp)))

    print("create output directories")
    # create output directories
    genelistsDir <<- glue("analysis/deseq2/{fct}/genelists")
    plotsDir <<- glue("analysis/deseq2/{fct}/plots")
    outDir <- list(genes=genelistsDir, plots=plotsDir)

    lapply(outDir, check_dir)

    print(glue("Output directories: {genelistsDir}"))
    print(glue("Output directories: {plotsDir}"))

    print("generating dds object")
    # Construct DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                colData = as.matrix(sampleinfo),
                                design = formula(glue("~ {fct}")))

    dds

    print("running deseq2")
    # Run DESeq2
    dds <- DESeq(dds)

    print("extracting comparisons results")
    # Extract results
    deseqDF <- lapply(cmp, get_deseq_df, dds, fct)
    deseqDF <- bind_cols(deseqDF)

    class(deseqDF)

    print("generating filtered gene lists")
    ## Export filtered output
    deg_list <- filterDEGs(degDF = deseqDF)
    saveRDS(deg_list, file = glue("{genelistsDir}/deseq2_results.RDS"))

    ## Exporting data

    print("exporting list")
    print(deg_list$Summary)
    # export DEG lists
    list_export(deg_list, annot, genelistsDir)

    print("exporting tables and plots")
    # export summary table
    write.table(deg_list$Summary,
        file = glue("{genelistsDir}/deg_summary.csv"),
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE)

    # export summary plot
    SavePlots(plots = deg_list[["Plot"]], 
        name = glue("{plotsDir}/DEG_summary"))

    # export upset plots (if there are more than one comparisons)
    upset_plots(deg_list, plotsDir)
}

## Set up factors to analyze
fct_comp <- names(sampleinfo)

lapply(fct_comp, deseq_per_factor)