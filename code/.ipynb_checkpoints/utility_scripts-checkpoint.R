#!/usr/bin/env Rscript

#### Helper script

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


#### Utility script for DEG analysis ####

generate_volcano <- function(df, l2fc="1", pval="0.05", topgenes="30"){
  
  df <- df %>% 
    mutate(diff_expressed = case_when(
        log2c > l2fc & padj < pval ~ "UP",
        log2c < -1*l2fc & padj < pval ~ "DOWN",
        TRUE ~ "NS"
        ))
  
  df$deg_label <- ifelse(df$gene_symbol %in% head(df[order(df$padj), "gene_symbol"], topgenes), df$gene_symbol, NA)
  
  vplot <- df %>%
    ggplot(aes(x = log2fc, y = -log10(pval), col = diff_expressed, label = deg_label)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 2) + 
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                        labels = c("Down-regulated", "Not significant", "Up-regulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Severe', #legend_title, 
        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
    geom_text_repel(max.overlaps = Inf) # To show all labels 

  vplot_obj <- list("df" = df, "plot" = vplot)

  return(vplot_obj)
}

lib_size_plot <- function() {


}

# Function to obtain all deseq outputs for each comparison
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

# Filter DEGs from deseq output
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

    data <- data[rownames(data) %in% unlist(DEGlistUPorDOWN[cmp]),]
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

# Export gene lists from deg_list output
list_export <- function(deg_list, annotations, genelistsDir){

  print(glue("Printing to outDir: {genelistsDir}"))
  deg_set <- c("Up", "Down")

  for(deg in deg_set){
    for(cmp in names(deg_list[[deg]])){

      genelist <- deg_list[[deg]][[cmp]] %>% as_tibble()
      colnames(genelist) <- "Gene ID"

      if(!is.null(deg_list[[deg]][[cmp]])){

        genelist <- genelist %>%
          left_join(annotations)
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

# Function to run GO analysis over each ontology
GOenrich <- function(deg_list, go_annot){

    results <- list()

    for(ont in names(go_annot)){

        go_ont <- go_annot[[ont]]
        go2gene <- go_ont[, c("GO ID", "Gene ID")]
        go2term <- go_ont[, c("GO ID", "GO_term")]

        x <- enricher(deg_list,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH", # holm, hochberg, bonferroni, BH, BY, fdr, none
            TERM2GENE = go2gene,
            TERM2NAME = go2term)

        results <- append(results, setNames(list(x), glue("{ont}")))

    }

    return(results)
}

# Function to export GO results
GO_results_export <- function(deg, allresults){

    for(cmp in names(allresults[[deg]])){

        mylist <- allresults[[deg]]

        merge <- data.frame(ID = character(),
                            Description = character(),
                            GeneRatio = character(),
                            BgRatio = character(),
                            pvalue = double(),
                            p.adjust = double(),
                            qvalue = double(),
                            geneID = character(),
                            Count = integer(),
                            stringsAsFactors = FALSE)

        for(ont in names(mylist[[cmp]])){

            if(!is.null(mylist[[cmp]][[ont]])){

            res <- mylist[[cmp]][[ont]]@result
            res$ontology <- ont

            merge <- dplyr::bind_rows(merge, res)
            rownames(merge) <- NULL
            }
        }

        write_csv(merge, 
                file = glue("{GODir}/{cmp}.{deg}.go.txt"),
                quote = "needed",
                col_names = TRUE
                )

    }
}

# Plotting KEGG pathways using Pathview
plot_kegg <- function(pathway_id, out_dir, fct, cmp){

  pathview::pathview(gene.data = ids,
    pathway.id = pathway_id,
    species = "btab",
    kegg.dir = out_dir,
    out.suffix = glue("{fct}_{cmp}"))

  # return(tryCatch(is.null(pview), error=function(e) NULL))

  system2("mv", args = c(glue("btab*.{fct}_{cmp}*"), KEGGDir_cmp))

}
