rm(list = ls())
#no sv*sex
#source("C:/Users/amadh/Downloads/rna_seq_functions(1).R")

library(edgeR)
library(org.Hs.eg.db)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(biomaRt)
library(readr)
library(ggplot2)
library(ggrepel)
library(sva)
library(DGEobj.utils)

setwd(paste0("C:/Users/amadh/Documents/RNA_SEQ/AMPAD"))

rosmap_unnorm <- read.delim(file ="AMP-AD_MSBB_MSSM_BM_10.raw_counts.tsv" )
# rnaseq_meta <- read.csv(file ="MSBB_assay_rnaSeq_metadata_corrected(1).csv")
pheno_data <- read.csv("MSBB_individual_metadata(1).csv")
rnaseq_cov <- read.csv("MSBB_RNAseq_covariates_November2018Update.csv")



rnaseq_cov <- rnaseq_cov[rnaseq_cov$BrodmannArea == "BM10", ]
rnaseq_cov <- rnaseq_cov[rnaseq_cov$Action == "OKay",]
IDS_df   <- rnaseq_cov[,c("sampleIdentifier", "individualIdentifier")]
rownames(rosmap_unnorm) <- rosmap_unnorm$Ensembl.ID
row_names <- rosmap_unnorm$Ensembl.ID
pheno_data <- pheno_data[pheno_data$individualID %in% IDS_df$individualIdentifier,]
IDS_df$dups <- duplicated(IDS_df$individualIdentifier)
IDS_df      <- IDS_df[IDS_df$dups == FALSE,]

IDS_df_chk <- IDS_df[IDS_df$sampleIdentifier %in% colnames(rosmap_unnorm),]
rosmap_unnorm  <- rosmap_unnorm[,colnames(rosmap_unnorm) %in% IDS_df_chk$sampleIdentifier]

rnaseq_cov     <- rnaseq_cov[rnaseq_cov$individualIdentifier %in% IDS_df_chk$individualIdentifier,]
rnaseq_cov$dup <- duplicated(rnaseq_cov$individualIdentifier)
rnaseq_cov     <- rnaseq_cov[rnaseq_cov$dup == FALSE,]
pheno_data     <- pheno_data[pheno_data$CERAD !=2,]
pheno_data$case  <- ifelse(pheno_data$CERAD >2, "AD", "Normal")

rnaseq_cov <- rnaseq_cov[rnaseq_cov$Action == "OKay",]
IDS_df   <- rnaseq_cov[,c("sampleIdentifier", "individualIdentifier")]

pheno_data <- pheno_data[pheno_data$individualID %in% IDS_df$individualIdentifier,]
rosmap_unnorm <- rosmap_unnorm[,2:ncol(rosmap_unnorm)]
rosmap_unnorm <- as.data.frame(rosmap_unnorm)

IDS_pheno <- IDS_df[IDS_df$individualIdentifier %in% pheno_data$individualID,]

counts_double <- apply(rosmap_unnorm,2,as.double)


counts_cpm <- convertCounts(counts_double, unit = "CPM")

rownames(counts_cpm) <- rownames(rosmap_unnorm)
counts_match_cpm <- counts_cpm[,colnames(counts_cpm) %in% IDS_pheno$sampleIdentifier]

keep <- rowSums(counts_match_cpm >= 10) > min(table(pheno_data$case, pheno_data$sex))
rosmap_final_count <- rosmap_unnorm[keep,]
rosmap_final_count <- rosmap_final_count[,colnames(rosmap_final_count) %in% IDS_pheno$sampleIdentifier]
IDS_pheno <- IDS_pheno[IDS_pheno$sampleIdentifier %in% colnames(rosmap_final_count),]

rownames(pheno_data) <- pheno_data$individualID
rownames(IDS_pheno) <- IDS_pheno$individualIdentifier

pheno_data <- merge(pheno_data,IDS_pheno, by = "row.names")
rownames(pheno_data) <- pheno_data$sampleIdentifier
id_order <- match(rownames(pheno_data), colnames(rosmap_final_count))
rosmap_final_count <- rosmap_final_count[,id_order]
pheno_data$age <- parse_number(pheno_data$ageDeath)
sva_df <- merge(t(rosmap_final_count), pheno_data[,c("case","sex","age")], by = "row.names")
rownames(sva_df) <- sva_df$Row.names
sva_df <- sva_df[,2:ncol(sva_df)]
colnames(sva_df) <- c(rownames(rosmap_final_count), "case","sex","age")
