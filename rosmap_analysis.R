setwd(paste0("C:/Users/axv851/Documents/RNA_SEQ/rosmap_file"))

library(edgeR)
library(org.Hs.eg.db)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(sva)
library(readr)
library(ggplot2)
library(cluster)
library(factoextra)
library(DGEobj.utils)

rosmap_counts <- read.delim("~/RNA_seq/rosmap_file/ROSMAP_all_counts_matrix (1).txt.gz", header=FALSE)
rosmap_meta <- read_csv("RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
rosmap_clinical <- read_csv("ROSMAP_clinical (1).csv")
rosmap_pheno <- read.delim("~/RNA_seq/rosmap_file/ROSMAP_Covariates_ages_censored.tsv", header=FALSE)

colnames(rosmap_pheno) <- rosmap_pheno[1,]
rosmap_pheno  <- rosmap_pheno[2:nrow(rosmap_pheno),]
colnames(rosmap_counts) <- rosmap_counts[1,]
rownames(rosmap_counts) <- rosmap_counts[,1]
rosmap_counts <- rosmap_counts[6:nrow(rosmap_counts),]
rosmap_counts <- rosmap_counts[,2:ncol(rosmap_counts)]


rosmap_pheno <- rosmap_pheno[rosmap_pheno$specimenID %in% colnames(rosmap_counts),]
rosmap_pheno$dup <- duplicated(rosmap_pheno$individualID)
rosmap_pheno     <- rosmap_pheno[rosmap_pheno$diagnosis %in% c("AD","CT"),]
rosmap_pheno$case <- ifelse(rosmap_pheno$diagnosis == "AD","AD","Normal")
rosmap_pheno$age <- parse_number(rosmap_pheno$age_death)

rosmap_counts <- rosmap_counts[,colnames(rosmap_counts) %in% rosmap_pheno$specimenID]
row_names <- rownames(rosmap_counts)
rosmap_counts <- apply(rosmap_counts,2,as.double)
rownames(rosmap_counts) <- row_names
rosmap_cpm    <- convertCounts(rosmap_counts, unit = "cpm")
keep <- rowSums(rosmap_cpm >= 1 ) >= min(table(rosmap_pheno$case, rosmap_pheno$sex))
rosmap_counts <- rosmap_counts[keep,]

rownames(rosmap_pheno) <- rosmap_pheno$specimenID
sva_df <- merge(t(rosmap_counts), rosmap_pheno[,c("case","sex","age")], by = "row.names")
rownames(sva_df) <- sva_df$Row.names
sva_df  <- sva_df[,2:ncol(sva_df)]