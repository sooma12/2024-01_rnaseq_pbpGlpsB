# deseq_output_analysis-2024-02-16.R

# BiocManager::install("DESeq2")

# Load libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)

data_wd="/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/rnaSeq/2024-01_rnaseq_pbpGlpsB/data"
setwd(data_wd)

DES_pbpG_df = read.csv("./DESeq_outputs/DES_pbpG_2024-02-16.csv", row.names = 1)

DES_lpsB_df = read.csv("./DESeq_outputs/DES_lpsB_2024-02-16.csv", row.names = 1)

#####
### Add functional annotations to DESeq output, so that we can select DEGs and pull KEGG/GO annotations
#####

# Read in functional annotations... `row.names = 2` sets row names to ACX60 locus tags
functional_annot_df <- read.csv("/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/REFERENCE/17978-mff_functional_annotations.csv", row.names = 2)

# functional annotations csv has some extraneous columns... grab only what we want.
selected_functional_annos <- functional_annot_df %>% dplyr::select(c("Global.Gene", "Global.Product", "Protein.ID", "Protein.names", "K.Number", "Tag1", "Tag2", "Tag3", "Category1", "Category2", "Category3", "Gene.ontology..GO.", "Gene.ontology..cellular.component.", "Gene.ontology..biological.process.","Gene.ontology..molecular.function."))

# Merge selected annotations into DESeq outputs
pbpG_DES_with_fxnannos <- left_join(as_tibble(DES_pbpG_df, rownames = "rn"), as_tibble(selected_functional_annos, rownames = "rn"), by = "rn")
lpsB_DES_with_fxnannos <- left_join(as_tibble(DES_lpsB_df, rownames = "rn"), as_tibble(selected_functional_annos, rownames = "rn"), by = "rn")

# Write out!
write.table(pbpG_DES_with_fxnannos, file = "./DESeq_processed/DES_d-pbpG_with_fxns_2024-02-26.csv", sep = ',', row.names = FALSE)
write.table(lpsB_DES_with_fxnannos, file = "./DESeq_processed/DES_d-lpsB_with_fxns_2024-02-26.csv", sep = ',', row.names = FALSE)

#just_names_and_l2fc_pbpG <- pbpG_DES_with_fxnannos %>% dplyr::select(c("rn", "log2FoldChange"))
#just_names_and_l2fc_lpsB <- lpsB_DES_with_fxnannos %>% dplyr::select(c("rn", "log2FoldChange"))
#write.table(just_names_and_l2fc_pbpG, file="./DESeq_processed/d-pbpG_loci_and_l2fc.tsv", sep = '\t', row.names = FALSE)
#write.table(just_names_and_l2fc_lpsB, file="./DESeq_processed/d-lpsB_loci_and_l2fc.tsv", sep = '\t', row.names = FALSE)
just_names_pbpG <- pbpG_DES_with_fxnannos %>% dplyr::select(c("rn"))
write.table(just_names_pbpG, file="./DESeq_processed/d-pbpG_loci.tsv", sep = '\t', row.names = FALSE)


#####
### Look at top upregulated and downregulated genes
#####
DES_pbpG_df_omitNA <- DES_pbpG_df %>% na.omit()
pbpG_downregulated <- head(DES_pbpG_df_omitNA[order(DES_pbpG_df_omitNA$log2FoldChange),], n=10L)
pbpG_upregulated <- tail(DES_pbpG_df_omitNA[order(DES_pbpG_df_omitNA$log2FoldChange),], n=10L)

functional_annot <- read.csv("/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/REFERENCE/17978-mff_functional_annotations.csv")
locus_to_gene <- functional_annot %>% select(c("Gene", "Global.Gene"))
row.names(locus_to_gene) <- locus_to_gene[,1]
locus_to_gene[,1] <- NULL

pbpG_downregulated_genes <- merge(pbpG_downregulated, locus_to_gene, by=0) %>% arrange(padj)
pbpG_upregulated_genes <- merge(pbpG_upregulated, locus_to_gene, by=0) %>% arrange(padj)

#####
### PlotMA
#####
# PlotMA needs: 3 columns (mean reads [numeric], logFC [numeric], and a logical vector for significance)

## pbpG
pbpG_to_MA_df <- as.data.frame(DES_pbpG_df) %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
pbpG_to_MA_df <- pbpG_to_MA_df %>% select(c("baseMean", "log2FoldChange", "padj"))
plotMA(pbpG_to_MA_df, ylim=c(-2, 2))

## lpsB TODO!

#####
### Enrichment analysis
#####

library(clusterProfiler)
library(KEGG.db)
functional_annot_df <- read.csv("/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/REFERENCE/17978-mff_functional_annotations.csv")


# Get list to enrich
sig_pbpG_degs <- DES_pbpG_df %>% filter(padj < 0.1 & abs(log2FoldChange) > 1)
sig_lpsB_degs <- DES_lpsB_df %>% filter(padj < 0.1 & abs(log2FoldChange) > 1)

annotations_for_pbpG_DEGs <- functional_annot_df %>% filter(functional_annot_df$Gene %in% row.names(sig_pbpG_degs))
annotations_for_lpsB_DEGs <- functional_annot_df %>% filter(functional_annot_df$Gene %in% row.names(sig_lpsB_degs))
common_DEGs <- inner_join(annotations_for_lpsB_DEGs, annotations_for_pbpG_DEGs)

# Convert gene IDs from genes list to K# ('K Number') from 17978-mff_functional_annotations.csv
genes_ko <- functional_annot_df$K.Number[match(row.names(sig_pbpG_degs), functional_annot_df$Gene)]
genes_ko <- genes_ko %>% na_if("") %>% na.omit()  # Remove any NAs or ""s from missing K#s

all_ab_kos <- locus_annotations$K.Number %>% na_if("") %>% na.omit()

ko_enrich <- enrichKEGG(gene = genes_ko, 
                        organism = "ko", 
                        keyType= "kegg",
                        universe = as.character(all_ab_kos), 
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1)

ko_enrich@result
GeneRatio <- ko_enrich@result[['GeneRatio']]
GeneRatio <- unname(sapply(GeneRatio, function(x) eval(parse(text = x))))
Description <- ko_enrich@result[['Description']]
Count <- ko_enrich@result[['Count']]
pValue <- ko_enrich@result[['pvalue']]
qValue <- ko_enrich@result[['qvalue']]
ko_result_tograph <- data.frame(GeneRatio, Description, Count, pValue, qValue)

select()

cat(row.names(sig_lpsB_degs), sep = '\n')
cat(row.names(sig_pbpG_degs), sep = '\n')
