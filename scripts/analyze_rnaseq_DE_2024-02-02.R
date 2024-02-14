# analyze_rnaseq_DE_2024-02-02.R
# DESeq2

# Installed these packages using Packrat:
# install.packages("tidyverse")
# install.packages("vsn") ## not for this R version
# install.packages("ggrepel")
# install.packages("DEGreport") ## not for this R version
# install.packages("pheatmap")
# BiocManager::install(pattern="apeglm")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
#.libpaths()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")
library(apeglm)
library(packrat)

# Set up packrat (package manager)
pr_dir <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/software'
setwd(pr_dir)
packrat::init()


# Load libraries
library(DESeq2)
library(ggplot2)
#library(vsn)
library(dplyr)
library(tidyverse)
library(ggrepel)
# library(DEGreport)
library(pheatmap)

wd <- "/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/"
setwd(wd)

## Read in raw data and prepare feature count table
feature_count <- read.table("./featurecounts/counts_noposttrim.txt", header=TRUE, row.names = 1)
# Grab the count data
data <- feature_count[,6:14]
# view(data)

## Set up column names in feature count data AND prep metadata table with strains and conditions
# Trim column names down to just the sample IDs
column_names <- colnames(data)
column_names <- sub("X.work.geisingerlab.Mark.rnaSeq.2024.01_rnaseq_pbpGlpsB.data.mapped.", "", column_names)
column_names <- sub("Aligned.sortedByCoord.out.bam", "", column_names)
column_names
colnames(data) <- column_names
# Use regex to get condition (mutation) from strain IDs
# ".+?(?=_)":   .+?  means "match any character (.) any number of times (+?)"
# (?=_): a positive lookahead: find where the character _ is matched
conditions <- str_extract(column_names, ".+?(?=_)")
# Use column_names and conditions to make metadata table
meta <- data.frame(column_names, conditions)
colnames(meta) <- c('id', 'condition')
meta  # Verify that IDs and conditions are as expected
# Write out metadata table (save!!)
#meta_filepath <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data'
#meta_file <- file.path(meta_filepath, 'metadata.txt')
#write.table(meta, meta_file, row.names=FALSE)  # TODO: Come back to this. Ok to have "id" label, or needs to be blank as in examples?  May want to change colnames before this.

## Load data, pre-filter low-count genes, and relevel to set WT as the reference
des_data <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
smallestGroupSize <- 3  # should be size of smallest group; I did 3 replicates
keep <- rowSums(counts(des_data) >= 10) >= smallestGroupSize  # keep data where count is >10 in all 3 samples
# relevel dds$condition to set WT as the reference
des_data$condition <- relevel(des_data$condition, ref = "WT")

dds <- DESeq(des_data)

# LFC shrinkage helps visualize and rank genes...
resultsNames(dds)  # get names of coefficients to shrink
## "Intercept" "condition_DlpsB_vs_WT" "condition_DpbpG_vs_WT"
# Use default shrinkage, which is apeglm method (Zhu, Ibrahim, and Love 2018)
resLFC_pbpG <- lfcShrink(dds, coef="condition_DpbpG_vs_WT", type="apeglm")
resLFC_lpsB <- lfcShrink(dds, coef="condition_DlpsB_vs_WT", type="apeglm")

# Count significant hits
sum(resLFC_pbpG$padj < 0.1, na.rm=TRUE)
sum(resLFC_lpsB$padj < 0.1, na.rm=TRUE)

# Check largest foldchanges from significant hits
resSigPbpG <- resLFC_pbpG[which(resLFC_pbpG$padj < 0.1), ]
#downregulated
head(resSigPbpG[order(resSigPbpG$log2FoldChange),])  # Verify pbpG (RS16895) is top hit
#upregulated
tail(resSigPbpG[order(resSigPbpG$log2FoldChange),])

resSigLpsB <- resLFC_lpsB[which(resLFC_lpsB$padj < 0.1), ]
#downregulated
head(resSigLpsB[order(resSigLpsB$log2FoldChange),])  # Verify lpsB (RS15950) is top hit
#upregulated
tail(resSigLpsB[order(resSigLpsB$log2FoldChange),])  # Verify lpsB (RS15950) is top hit

# MA plot
plotMA(resLFC_pbpG, ylim=c(-1, 1))
plotMA(resLFC_lpsB, ylim=c(-1, 1))

# Check PCA plot
#rld <- rlog(des, blind=TRUE)
#plotPCA(rld, intgroup="condition") + geom_text(aes(label=name))

# Plot the gene with the smallest value 
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# plotCounts(dds, gene="gene-ACX60_RS15100", intgroup="condition")

# Save data
#file_name = 'NAME/OF/FILE'
#write.table(<results/go/here>, file_name, quote=FALSE)
# Can also pre-filter for significant hits
#padj.cutoff <- 0.05 # False Discovery Rate cutoff
#significant_results <- res[which(res$padj < padj.cutoff),]
#write.table(significant_results, file_name, quote=FALSE)

#get session information, including package versions... can save to a variable, then write out.
#print(sessionInfo())


# Save current analyses
pbpGresOrdered <- resLFC_pbpG[order(resLFC_pbpG$pvalue),]
write.csv(as.data.frame(pbpGresOrdered), file="pbpG_results_2024-02-07.csv")
lpsBresOrdered <- resLFC_lpsB[order(resLFC_lpsB$pvalue),]
write.csv(as.data.frame(lpsBresOrdered), file="lpsB_results_2024-02-07.csv")
