library(edgeR)
library(org.Hs.eg.db)
require(DESeq2)
library(pheatmap)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(readr)
library(sva)
library(DGEobj.utils)
Mayo_gene_all_counts_matrix_clean <- read.delim("C:/Users/amadh/Downloads/Mayo_gene_all_counts_matrix_clean.txt", header=FALSE)
RNAseq_Harmonization_Mayo_combined_metadata <- read_csv("C:/Users/amadh/Downloads/RNAseq_Harmonization_Mayo_combined_metadata.csv")

RNAseq_Harmonization_Mayo_combined_metadata <- RNAseq_Harmonization_Mayo_combined_metadata[RNAseq_Harmonization_Mayo_combined_metadata$diagnosis %in% c("Alzheimer Disease", "control"),]
table(RNAseq_Harmonization_Mayo_combined_metadata$tissue)
RNAseq_Harmonization_Mayo_combined_metadata <- RNAseq_Harmonization_Mayo_combined_metadata[RNAseq_Harmonization_Mayo_combined_metadata$tissue == "cerebellum",]
colnames(Mayo_gene_all_counts_matrix_clean) <- Mayo_gene_all_counts_matrix_clean[1,]
rownames(Mayo_gene_all_counts_matrix_clean) <- Mayo_gene_all_counts_matrix_clean[,1]
Mayo_gene_all_counts_matrix_clean <- Mayo_gene_all_counts_matrix_clean[,2:ncol(Mayo_gene_all_counts_matrix_clean)]
Mayo_gene_all_counts_matrix_clean <- Mayo_gene_all_counts_matrix_clean[2:nrow(Mayo_gene_all_counts_matrix_clean),]

Mayo_gene_all_counts_matrix_clean <- Mayo_gene_all_counts_matrix_clean[5:nrow(Mayo_gene_all_counts_matrix_clean),]

mayo_gene_cer <- Mayo_gene_all_counts_matrix_clean[,colnames(Mayo_gene_all_counts_matrix_clean) %in% RNAseq_Harmonization_Mayo_combined_metadata$specimenID]
pheno_data <- RNAseq_Harmonization_Mayo_combined_metadata
pheno_data$case <- ifelse(pheno_data$diagnosis == "control", "Normal", "AD")
pheno_data$case <- as.factor(pheno_data$case)
pheno_data$sex <- as.factor(pheno_data$sex)
pheno_data$case <- relevel(pheno_data$case, ref = "Normal")
pheno_data <- as.data.frame(pheno_data)
rownames(pheno_data) <- pheno_data$specimenID



pheno_data <- pheno_data[rownames(pheno_data) %in% colnames(mayo_gene_cer),]
pheno_data <- pheno_data[order(as.character(row.names(pheno_data))),]
mayo_gene_cer <- mayo_gene_cer[,order(as.character(colnames(mayo_gene_cer)))]
order_check <- rownames(pheno_data) == colnames(mayo_gene_cer)
table(order_check)
mayo_gene_cer <- apply(mayo_gene_cer,2,as.double)
rownames(mayo_gene_cer) <- rownames(Mayo_gene_all_counts_matrix_clean)
mayo_gene_cpm <- convertCounts(mayo_gene_cer, unit = "cpm")
keep <- rowSums(mayo_gene_cpm >= 1) > min(table(pheno_data$case, pheno_data$sex))
mayo_gene_cer <- mayo_gene_cer[keep,]
pheno_data$age <- parse_number(pheno_data$ageDeath)

sva_df <- merge(t(mayo_gene_cer), pheno_data[,c("case","sex","age")], by = "row.names")
rownames(sva_df) <- sva_df$Row.names
sva_df <- sva_df[,2:ncol(sva_df)]

