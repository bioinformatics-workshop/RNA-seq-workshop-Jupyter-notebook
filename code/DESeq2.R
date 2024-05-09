#!/usr/bin/env Rscript

# Load libaries
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(EnhancedVolcano)
library(ggfortify)
library(glue)
library(plotly)

# Load count matrix file from featurecounts and clean up the header
counts <- read_tsv("analysis/featurecounts/featurecounts_full.txt",
          col_names=TRUE, skip=1)

new_name <- str_extract(names(counts)[endsWith(names(counts), "bam")],
            "(?<=star/).*(?=_Aligned)")

names(counts)[endsWith(names(counts), "bam")] <- new_name

counts <- counts %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

# Load sample info and metadata

sampleinfo <- read_csv("raw/metadata.csv", col_names=TRUE)

sampleinfo <- sampleinfo %>%
  mutate(genotype = factor(genotype, levels = c("wt", "miR163_mut"))) %>%
  arrange(samplename) %>%
  column_to_rownames(var = "samplename") %>%
  as.matrix()

write.table(as.data.frame(counts),
  file = "analysis/deseq2/genelists/count_matrix.csv",
  sep = ",",
  quote = FALSE,
  row.names = F)

# Load annotation file
annot <- read_csv("genome/ath_gene_annot.csv", col_names=TRUE)

# Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleinfo,
                              design = ~ genotype)

# Reorder factor level
dds$genotype<- relevel(dds$genotype, ref = "wt")

# Pre-filter low counts
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)

normalized.counts <- as.data.frame(counts(dds, normalized=TRUE)) %>%
  as_tibble(rownames = "locus") %>%
  left_join(., annot, c('locus'='tair_locus'))

# write to file
write.table(normalized.counts,"analysis/deseq2/genelists/DESeq2_normalizedcounts.txt", 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE)

############## Pre-Processing Quality Assessment (START) ##############

## Library sizes bar plot

librarySizes <- colSums(counts)

pdf("analysis/deseq2/plots/lib-size.barplots.pdf", width = 7, height = 7)
barplot(librarySizes,
        names=names(librarySizes),
        las=2,
        main="Barplot of library sizes")
abline(h=20e6, lty=2)
dev.off()

png("analysis/deseq2/plots/lib-size.barplots.png")
barplot(librarySizes,
        names=names(librarySizes),
        las=2,
        main="Barplot of library sizes")
abline(h=20e6, lty=2)
dev.off()

## Count distribution Boxplots

logcounts <- log2(counts + 1) # no normalization
title_boxplot <- paste("Count Distribution Boxplots")

pdf("analysis/deseq2/plots/count-dist.boxplots.raw.pdf", width = 7, height = 7)
boxplot(logcounts,
        main=title_boxplot,
        xlab="",
        ylab="Log2(Counts+1)",
        las=2,
        col=c(rep("red",3), rep("blue",3)))
abline(h=median(as.matrix(logcounts)), col="blue")
dev.off()

png("analysis/deseq2/plots/count-dist.boxplots.raw.png")
boxplot(logcounts,
        main=title_boxplot,
        xlab="",
        ylab="Log2(Counts+1)",
        las=2,
        col=c(rep("red",3),rep("blue",3)))
abline(h=median(as.matrix(logcounts)), col="blue")
dev.off()

rlogcounts <- rlog(counts) #normalization
title_boxplot_rlog <- paste("Count Distribution Boxplots (rlog)")

pdf("analysis/deseq2/plots/count-dist.boxplots.rlog.pdf", width = 7, height = 7)
boxplot(rlogcounts,
        main=title_boxplot_rlog,
        xlab="",
        ylab="rlog(Counts)",
        las=2,
        col=c(rep("red",3),rep("blue",3)))
dev.off()

png("analysis/deseq2/plots/count-dist.boxplots.rlog.png")
boxplot(rlogcounts,
        main=title_boxplot_rlog,
        xlab="",
        ylab="rlog(Counts)",
        las=2,
        col=c(rep("red",3),rep("blue",3)))
dev.off()


## Correlation heatmaps

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
rld <- rlog(dds, blind=TRUE)
corr_samps <- cor(as.matrix(assay(rld)))

pdf("analysis/deseq2/plots/cor-heatmap.pdf", width = 7, height = 7)
heatmap.2(corr_samps,
          trace = "none",
          col = colors,
          main = "Sample correlations")
dev.off()

png("analysis/deseq2/plots/cor-heatmap.png")
heatmap.2(corr_samps,
          trace = "none",
          col = colors,
          main = "Sample correlations")
dev.off()

## Principal Component Analysis

# autoplot requires the ggfortify library

pcDat <- prcomp(t(rlogcounts))
title_pca <- paste("PC1 vs PC2, All Samples")

pdf("analysis/deseq2/plots/PCA.pdf", width = 7, height = 7)
autoplot(pcDat,
         main=title_pca,
         data = sampleinfo,
         colour="genotype",
         shape=FALSE)
dev.off()

png("analysis/deseq2/plots/PCA.png")
autoplot(pcDat,
         main=title_pca,
         data = sampleinfo,
         colour="genotype",
         shape=FALSE)
dev.off()

## Plotly version of PCA plot
p <- autoplot(pcDat,
         main=title_pca,
         data = sampleinfo,
         colour = "genotype",
         label = F,
         label.repel = F,
         size = 8,
         group = "time",
         shape = "time")

plt <- plotly::ggplotly(p)
htmlwidgets::saveWidget(plt, "analysis/deseq2/plots/PCA.html")

############## Pre-Processing Quality Assessment (END) ##############


############## Differential Expression (START)##############

# Calculate DE for all contrasts
dds <- DESeq(dds)

# Extract different contrast (i.e., comparisons)
contrast_list <- resultsNames(dds)
contrast_list <- contrast_list[contrast_list != "Intercept"] # remove Intercept




# Function to extract info for each contrast

deseq2_deg <- function(contrast_name){
  # extracting DEG results
  res <- results(dds, name=contrast_name)
  res05 <- results(dds, name=contrast_name, lfcThreshold = 1, alpha = 0.05)
  res05$padj <- p.adjust(res05$pvalue, method="BH")
  res05 <- res05[!is.na(res05$padj),]

  res_filt <- res[ which(res$padj < 0.1), ] # total DEGs
  res05_filt <- res05[ which(res05$padj < 0.05), ] # total DEGs
  up <- subset(res_filt, log2FoldChange > 0) # upDEGs (LFC > 0, p<0.1)
  dn <- subset(res_filt, log2FoldChange < 0) # dnDEGs (LFC < 0, p<0.1)
  up_05 <- subset(res05_filt, log2FoldChange > 1) # upDEGs (LFC > 1, p<0.05)
  dn_05 <- subset(res05_filt, log2FoldChange < -1) # dnDEGs (LFC< -1, p<0.05)

  # generates summary file
  summary_file <- c(contrast_name,
                    nrow(res_filt),
                    nrow(up),
                    nrow(dn),
                    nrow(res05_filt),
                    nrow(up_05),
                    nrow(dn_05))

  # export lists
  degs <- subset(res, padj < 0.1) %>%
          as_tibble(rownames = "locus") %>%
          left_join(., annot, c('locus'='tair_locus'))

  degs_05 <- subset(res05, padj < 0.05, log2FoldChange > 1) %>%
             as_tibble(rownames = "locus") %>%
             left_join(., annot, c('locus'='tair_locus'))

  write.table(degs,
              file = paste0("analysis/deseq2/genelists/",contrast_name,"_DEGs.LFC0_p1.txt"),
              quote = FALSE,
              sep= "\t",
              row.names = FALSE)
  write.table(degs_05,
              file = paste0("analysis/deseq2/genelists/",contrast_name,"_DEGs.LFC1_p05.txt"),
              quote = FALSE,
              sep= "\t",
              row.names = FALSE)

  # post-processing plots

  # volcano plots with EnhancedVolcano

  subtitle1 <- paste("LFC > 0, FDR < 0.1")
  subtitle2 <- paste("LFC > 1, FDR < 0.05")

  EnhancedVolcano(res,
      lab = rownames(res),
      x = 'log2FoldChange',
      y = 'pvalue',
      title = contrast_name,
      subtitle = subtitle1,
      pCutoff = 0.1,
      pCutoffCol = 'padj', # p-value
      FCcutoff = 0,
      pointSize = 4.0,
      labSize = 6.0,
      drawConnectors = TRUE,
      widthConnectors = 0.75)

  ggsave(paste0("analysis/deseq2/plots/",contrast_name,"-volcano-LFC0-FDR1.pdf"), width = 7, height = 7)
  ggsave(paste0("analysis/deseq2/plots/",contrast_name,"-volcano-LFC0-FDR1.png"), scale = 2)

  EnhancedVolcano(res05,
      lab = rownames(res05),
      x = 'log2FoldChange',
      y = 'pvalue',
      title = contrast_name,
      subtitle = subtitle2,
      pCutoff = 0.05,
      pCutoffCol = 'padj', # p-value
      FCcutoff = 1,
      pointSize = 4.0,
      labSize = 6.0,
      drawConnectors = TRUE,
      widthConnectors = 0.75)

  ggsave(paste0("analysis/deseq2/plots/",contrast_name,"-volcano-LFC1-FDR05.pdf"), width = 7, height = 7)
  ggsave(paste0("analysis/deseq2/plots/",contrast_name,"-volcano-LFC1-FDR05.png"), scale = 2)

  return(summary_file)
}

# Run function for all contrast
sum_file <- sapply(contrast_list, deseq2_deg)

sum_file_tbl <- sum_file %>%
  t() %>%
  as_tibble() %>%
  filter(V1 != "Intercept") %>%
  select("Comparison"=V1,
         "Total DEGs (p < 0.1)"=V2,
         "Up DEGs (p < 0.1, LFC > 0)"=V3,
         "Down DEGs (p < 0.1, LFC < 0)"=V4,
         "Total DEGs (p < 0.05)"=V5,
         "Up DEGs (p < 0.05, LFC > 1)"=V6,
         "Down DEGs (p < 0.05, LFC < -1)"=V7) %>%
  mutate(across(2:7, as.integer))

write.table(sum_file_tbl,
            file = "analysis/deseq2/genelists/DESeq2.DEG_summary.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

