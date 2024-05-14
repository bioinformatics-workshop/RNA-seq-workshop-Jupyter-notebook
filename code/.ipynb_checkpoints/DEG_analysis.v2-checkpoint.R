#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(DESeq2)
library(glue)
library(UpSetR)

# Load annotation file
# bta_to_loc <- read_tsv("code/bta_to_loc.annotation.txt")

# bta2loc_filtered <- bta_to_loc %>%
#   group_by(query_id) %>%
#   filter(eval == min(eval) & identity == max(identity)) %>%
#   ungroup()

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

# filter DEGs
filterDEGs <- function(degDF, FDR = 0.05, Fold = 2) {

  pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
  log2FC <- degDF[, grep("_log2FC$", colnames(degDF)), drop = FALSE]

  ## DEGs that are up or down regulated
  pf <- pval <= FDR & (log2FC >= log2(Fold) | log2FC <= -log2(Fold))
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop = FALSE]), simplify = FALSE)

  # all_plt <- upset(fromList(DEGlistUPorDOWN),
  #   order.by = "freq",
  #   nsets = length(names(DEGlistUPorDOWN)),
  #   mainbar.y.label = "Number of DEGs Overlapped",
  #   sets.x.label = "Number of DEGs",
  #   group.by = "sets",
  #   mb.ratio = c(0.7,0.3),
  # print(all_plt)

  ## DEGs that are up regulated
  pf <- pval <= FDR & log2FC >= log2(Fold)
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)

  ## DEGs that are down regulated
  pf <- pval <= FDR & log2FC <= -log2(Fold)
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)
  df <- data.frame(Comparisons=names(DEGlistUPorDOWN), Counts_Up_or_Down=sapply(DEGlistUPorDOWN, length), Counts_Up=sapply(DEGlistUP, length), Counts_Down=sapply(DEGlistDOWN, length))
#  resultlist <- list(UporDown=DEGlistUPorDOWN, Up=DEGlistUP, Down=DEGlistDOWN, Summary=df)
  
  # summary plots
    df_plot <- data.frame(Comparisons=rep(as.character(df$Comparisons), 2), Counts=c(df$Counts_Up, df$Counts_Down), Type=rep(c("Up", "Down"), each=length(df[,1])))
    plt <- df_plot %>%
      dplyr::mutate(Counts = case_when(
        Type == "Down" ~ -1*Counts,
        TRUE ~ Counts
      )) %>%
    ggplot(aes(x = Comparisons, y = Counts, fill = Type)) +
    geom_col() +
  #  geom_text(aes(label = abs(Counts)), hjust = 0.01) +
    geom_hline(yintercept = 0, color = "black") +
    labs(y = "Number of DEGs",
        x = "",
        title = glue("DEG Counts (Fold = {Fold} & FDR = {FDR})")) +
  #  coord_cartesian(ylim=c(ceiling(-ylim*1.2),ceiling(ylim*1.2))) +
    coord_flip() +
    facet_grid(~Type, scales = "free_x") +
    scale_y_continuous(expand = c(0,0),
                      labels = function(x) abs(x)) +
    # theme_minimal() +
    # ylim(-4000, 4000) +

    theme(legend.position = "none",
          panel.spacing.x = unit(0, "pt"),
          strip.background = element_rect(color = "black", fill = "grey"),
          strip.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14))
      print(plt)

    # upset plot Up DEG
    up_plt <- upset(fromList(DEGlistUP),
      order.by = "freq",
      nsets = length(names(DEGlistUP)),
      mainbar.y.label = "Number of DEGs Overlapped",
      sets.x.label = "Number of DEGs",
      group.by = "sets",
      mb.ratio = c(0.7,0.3),
      )
    print(up_plt)

    # upset plot Down DEG
    down_plt <- upset(fromList(DEGlistDOWN),
      order.by = "freq",
      nsets = length(names(DEGlistDOWN)),
      mainbar.y.label = "Number of DEGs Overlapped",
      sets.x.label = "Number of DEGs",
      group.by = "sets",
      mb.ratio = c(0.7,0.3),
      )
    print(down_plt)

    # upset plot All DEG
    all_plt <- upset(fromList(DEGlistUPorDOWN),
      order.by = "freq",
      nsets = length(names(DEGlistUPorDOWN)),
      mainbar.y.label = "Number of DEGs Overlapped",
      sets.x.label = "Number of DEGs",
      group.by = "sets",
      mb.ratio = c(0.7,0.3),
    print(all_plt)

  plt_list <- list(DEG_summary=plt,
                  upset_up=up_plt,
                  upset_down=down_plt,
                  upset_all=all_plt)

  resultlist <- list(UporDown=DEGlistUPorDOWN,
                      Up=DEGlistUP,
                      Down=DEGlistDOWN,
                      Summary=df,
                      Plots=plt_list)


  return(resultlist)

}

# Export gene lists
list_export <- function(deg_list){

  print(glue("Printing to outDir: {genelistsDir}"))
  deg_set <- c("Up", "Down")

  for(deg in deg_set){
    for(cmp in names(deg_list[[deg]])){

      genelist <- deg_list[[deg]][[cmp]] %>% as_tibble()
      colnames(genelist) <- "Gene ID"

      if(!is.null(deg_list[[deg]][[cmp]])){

        genelist <- genelist %>%
          left_join(annot)
      }

      write_csv(genelist, 
                  file = glue("{genelistsDir}/{cmp}.{deg}DEG.txt"),
                  quote = "needed",
                  col_names = TRUE
                  )
    }
  }
}


plots_export <- function(deg_list){

  for(plot in names(deg_list[["Plots"]])){

    SavePlots(plots = deg_list[["Plots"]][[plot]], 
        name = glue("{plotsDir}/{plot}"))

  }

}




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
  deseqDF <- lapply(cmp, get_deseq_df, dds = dds, fct = fct)
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
    list_export(deg_list, annotations = annot, genelistsDir = genelistsDir)

  print("exporting tables and plots")
  # export summary table
    write.table(deg_list$Summary,
              file = glue("{genelistsDir}/deg_summary.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)

  # export plot
  SavePlots(plots = deg_list[["Plot"]], 
      name = glue("{plotsDir}/DEG_plot"))

  # export upset plots
  upset_plots(deg_list, plotsDir)
}

## Set up factors to analyze
# fct_comp <- c("fct_infection",
#               "fct_genotype",
#               "fct_sex",
#               "fct_substrate",
#               "fct_infection_substrate",
#               "fct_infection_genotype",
#               "fct_infection_sex",
#               "fct_genotype_sex",
#               "fct_infection_genotype_sex",
#               "fct_infection_genotype_substrate",
#               "fct_infection_sex_substrate",
#               "fct_genotype_sex_substrate",
#               "fct_all")

fct_comp <- names(sampleinfo)

lapply(fct_comp, deseq_per_factor)

fct <- "fct_all"
