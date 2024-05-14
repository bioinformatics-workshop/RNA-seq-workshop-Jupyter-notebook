#!/usr/bin/env Rscript

## R functions to assist with the data processing


#### Helper function

# SavePlots - For fast saving of ggplots to PNG and PDF formats
SavePlots <- function(plots, name, width=13, height=8){   

  ggsave(plots, filename=paste0(name,".pdf"), width = width, height = height) 
  ggsave(plots, filename=paste0(name,".png"), width = width, height = height) 

}

# creates output directories
check_dir <- function(dir){

    if( !dir.exists(dir) ){
      dir.create(dir, recursive = TRUE)
    }
}


#### Utility functions for DEG analysis ####

# Merge featurecount files into one
fcount_merge <- function(dir = "analysis/featurecounts"){

    list_of_files <- list.files(dir,
      pattern = "\\.fcnts.txt$",
      full.names = TRUE)

    extract_data_from_file <- function(file){
      data <- read_delim(file,
              id = "sample",
              comment = "#",
              delim = "\t",
              col_names = TRUE) %>%
              dplyr::select(Geneid, Counts = last_col(), sample) %>%
              mutate(sample = str_extract(sample,"(?<=featurecounts/).*(?=.fcnts)"))
      
      return(data)
    }

    alldata <- map_dfr(list_of_files, extract_data_from_file)

    alldata <- alldata %>%
        pivot_wider(
        names_from = sample,
        values_from = Counts,
        names_sort = TRUE) %>%
        column_to_rownames("Geneid")

    print(alldata)
}


# Function to obtain all deseq2 outputs for each comparison
get_deseq_df <- function(cmp, dds, fct) {
  deseq_df <- data.frame(row.names = rownames(counts))
  res <- DESeq2::results(dds, contrast = c(fct, cmp))
  ## Set NAs to reasonable values to avoid errors in downstream filtering steps
  res[is.na(res[, "padj"]), "padj"] <- 1
  res[is.na(res[, "log2FoldChange"]), "log2FoldChange"] <- 0
  deg <- as.data.frame(res)
  colnames(deg)[colnames(deg) %in% c("log2FoldChange", "padj")] <- c("log2FC", "FDR")
  colnames(deg) <- paste(paste(cmp, collapse = "-"), colnames(deg), sep = "_")
  deseq_df <- cbind(deseq_df, deg[rownames(deseq_df), ])

  return(deseq_df)
}


# Filter DEGs from deseq2 output
filterDEGs <- function(degDF, FDR = 0.05, Fold = 2){

  pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
  log2FC <- degDF[, grep("_log2FC$", colnames(degDF)), drop = FALSE]

  ## DEGs that are up or down regulated
  pf <- pval <= FDR & (log2FC >= log2(Fold) | log2FC <= -log2(Fold))
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop = FALSE]), simplify = FALSE)

  extract_deseq_data <- function(cmp, filter = TRUE){

    data <- degDF[, grep(cmp, colnames(degDF)), drop = FALSE]

    if (filter == TRUE) {

    data <- data[rownames(data) %in% unlist(DEGlistUPorDOWN[cmp]), ]
    colnames(data) <- gsub(glue("{cmp}_"), "", colnames(data))

    } else {

    colnames(data) <- gsub(glue("{cmp}_"), "", colnames(data))

    }

    return(data)
    }

  ## DESeq2 results for DEGs only
  deg_data <- lapply(colnames(pf), extract_deseq_data)
  names(deg_data) <- colnames(pf)

  ## DESeq2 results for all genes
  all_data <- lapply(colnames(pf), extract_deseq_data, filter = FALSE)
  names(all_data) <- colnames(pf)

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

  ## Summary dataframe
  df <- data.frame(Comparisons=names(DEGlistUPorDOWN),
        Counts_Up_or_Down=sapply(DEGlistUPorDOWN, length), Counts_Up=sapply(DEGlistUP, length),
        Counts_Down=sapply(DEGlistDOWN, length))

  ## summary plots
  df_plot <- data.frame(Comparisons=rep(as.character(df$Comparisons), 2), Counts=c(df$Counts_Up, df$Counts_Down), Type=rep(c("Up", "Down"), each=length(df[,1])))
  
  plt <- df_plot %>%
    dplyr::mutate(Counts = case_when(
      Type == "Down" ~ -1*Counts,
      TRUE ~ Counts
    )) %>%
  ggplot(aes(x = Comparisons, y = Counts, fill = Type)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "black") +
  labs(y = "Number of DEGs",
      x = "",
      title = glue("DEG Counts (Fold = {Fold} & FDR = {FDR})")) +
  coord_flip() +
  facet_grid(~Type, scales = "free_x") +
  scale_y_continuous(expand = c(0,0),
                    labels = function(x) abs(x)) +
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "pt"),
        strip.background = element_rect(color = "black", fill = "grey"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

  print(plt)


    resultlist <- list(AllData = all_data,
                      DEGData = deg_data,
                      UporDown = DEGlistUPorDOWN,
                      Up = DEGlistUP,
                      Down = DEGlistDOWN, 
                      Summary = df,
                      Plot = plt)

    return(resultlist)

}

# Generate Enhanced Volcano plots
e_volcano_plot <- function(cmp, dds, fct){

  res <- DESeq2::results(dds, contrast = c(fct, cmp))

  cond1 <- cmp[[1]]
  cond2 <- cmp[[2]]

  title <- glue("Factor: {fct} \n Comparison: {cond1} vs {cond2}") 
  subtitle <- paste("LFC > 1, FDR < 0.05")
  
  plt <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = title,
    subtitle = subtitle,
    pCutoff = 0.05,
    pCutoffCol = 'padj', # p-value
    FCcutoff = 0,
    pointSize = 4.0,
    labSize = 6.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75)

  SavePlots(plt, name = glue("{plotsDir}/Enhancedvolcano_{cond1}-vs-{cond2}"))
}


# export for Ath annotation only
list_export <- function(deg_list, annotations, genelistsDir){

  print(glue("Printing to outDir: {genelistsDir}"))
  deg_set <- c("Up", "Down")

  for(deg in deg_set){
    for(cmp in names(deg_list[[deg]])){

      genelist <- deg_list[[deg]][[cmp]] %>% as_tibble()
      colnames(genelist) <- "Gene ID"

      if(!is.null(deg_list[[deg]][[cmp]])){

        genelist <- genelist %>%
          left_join(annotations, by = c("Gene ID" = "tair_locus"))
      }

      write_csv(genelist, 
                  file = glue("{genelistsDir}/{cmp}.{deg}DEG.annotated.csv"),
                  quote = "needed",
                  col_names = TRUE
                  )
    }
  }

  for(cmp in names(deg_list[["DEGData"]])){

    degdata <- deg_list[["DEGData"]][[cmp]] %>%
    as_tibble(rownames = "gene_id")

    write_csv(degdata,
              file = glue("{genelistsDir}/{cmp}.deseq2_output.csv"),
              quote = "needed",
              col_names = TRUE)

  }
  
}


# Export plots from deg_list object
upset_plots <- function(deg_list, plotsDir) {

  plot_list <- c("UporDown", "Up", "Down")

  for(degs in plot_list) {
    
    if(length(names(deg_list[[degs]])) > 1){
    my_plot <- UpSetR::upset(UpSetR::fromList(deg_list[[degs]]),
      order.by = "freq",
      nsets = length(names(deg_list[[degs]])),
      mainbar.y.label = "Number of DEGs Overlapped",
      sets.x.label = "Number of DEGs",
      group.by = "sets",
      mb.ratio = c(0.7,0.3)
      ) 

    pdf(file = glue("{plotsDir}/upset_{degs}.pdf"))
    print(my_plot)
    dev.off()

    png(file = glue("{plotsDir}/upset_{degs}.png"))
    print(my_plot)
    dev.off()
    }
  }  
}

# PCA plot
generate_pca <- function(dds, fct) {

  rld <- rlog(dds)

  pltdata <- plotPCA(rld, intgroup = fct, returnData = TRUE)
  percentVar <- round(100 * attr(pltdata, "percentVar"))
  plt <- ggplot(pltdata, aes(PC1, PC2, color=group, label=name)) +
          geom_point(size=3) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          coord_fixed()

  pltly <- plotly::ggplotly(plt)
  htmlwidgets::saveWidget(pltly, glue("{plotsDir}/PCA.html"))

  SavePlots(plt, glue("{plotsDir}/PCA"))

}
