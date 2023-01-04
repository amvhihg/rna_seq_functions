#' Calculate pca decomposition of an expression matrix and graph a scatter plot of pca values colored by phenotype of choice
#' 
#' @param expr_m an normalized and variance stabilized expression matrix
#' @param exp_var A character signifying the name of the phenotype used for coloring scatterplot points in the pca plot
#' 
#' @return list containing 1.) scatter plot of the first two principal components 2.) output from the functino prcomp on the expression matrix 3.)
#' mean and standard deviations of the first two principal components 4.) A dataframe containing the first two principal components and the phenotype of interest

pca_data_out <- function(expr_m, exp_var){
  

  n_cols            <- ncol(expr_m) - 3 # Final two columns contain the expression matrix
  
  gene_exp_mat  <- t(as.matrix(expr_m[,1:n_cols]))  # gene expression matrix
  gene_exp_mat  <- apply(gene_exp_mat,2,as.numeric) # convert the values to numeric 
  gene_exp_mat  <- t(gene_exp_mat)                  # samples as columns and features as rows  
  pca           <- prcomp( gene_exp_mat, center = TRUE, scale. = TRUE ) # build pca object
  pca_df        <- as.data.frame(cbind(pca$x[,1],pca$x[,2],expr_m[exp_var])) # create a dataframe with the first two principal components and the phenotype variable
  colnames(pca_df) <- c("PC1", "PC2","status") 
  pca_df$status <- as.factor(pca_df$status)
  pca_df$pat    <- rownames(pca_df)
  # create a ggplot object with the phenotype column as the color 
  hor1          <- sd(pca$x[,1])
  hor2          <- sd(pca$x[,2])
  ver1          <- mean(pca$x[,1])
  ver2          <- mean(pca$x[,2])
  exp_pca       <- ggplot(pca_df, aes(x = PC1,y = PC2, label = rownames(pca_df))) + geom_point(aes(color = status)) 
  exp_pca       <- exp_pca + geom_vline( xintercept = ver1 + 3*hor1 ) + geom_vline( xintercept = ver1 - 3*hor1)    + geom_text_repel()# add vertical lines to cover +- 3 standard deviations
  exp_pca       <- exp_pca + geom_hline( yintercept = ver2 + 3*hor2)  +  geom_hline(yintercept = ver2 - 3*hor2) # add horizontal lines to cover +-3 standard deviations

  
  
  return(list(exp_pca, pca,list(hor1,hor2,ver1,ver2),pca_df))
  
}
#' Function to estimate surrogate variable values from a rnaseq generated expression matrix
#' 
#' @param edat A dataframe with transpose of an expression matrix merged with case, sex and age
#' 
#' @return list composed of output from svaseq function from package "sva"
sva_cleaning_age <- function(edat){


  mod_data   <- model.matrix( ~1 + case + sex+ case *sex + age + age*sex, data = edat)
  mod0       <- model.matrix(~1  +sex + age + age *sex , data = edat)
  n_cols     <- ncol(edat) - 3
  p_cols     <- ncol(edat) - 1
  expr_ae_ad <- apply(edat[,1:n_cols],2,as.double)
  #expr_ae_ad <- as.matrix(eset_ae_ad_pb[,1:2063])
  #expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)
  
  n.survar   <- num.sv(t(expr_ae_ad),mod_data, method = "leek") 
  
  svobj  <- svaseq(t(expr_ae_ad), mod_data, mod0,numSVmethod = "leek")
  
  return(svobj)
}
sva_cleaning_mod <- function(edat){
  # edat is the expression matrix with sex and case as the last two columns 
  mod_data   <- model.matrix( ~1 + case + sex+ case *sex, data = edat)
  mod0       <- model.matrix(~1  +sex  , data = edat)
  n_cols     <- ncol(edat) - 3
  p_cols     <- ncol(edat) - 1
  expr_ae_ad <- apply(edat[,1:n_cols],2,as.double)
  #expr_ae_ad <- as.matrix(eset_ae_ad_pb[,1:2063])
  #expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)
  
  n.survar   <- num.sv(t(expr_ae_ad),mod_data, method = "leek") 
  
  svobj  <- svaseq(t(expr_ae_ad), mod_data, mod0,numSVmethod = "leek")
  
  return(svobj)
}

#' Map ensembl gene ids to hgnc gene symbols and chromosome 
#' 
#' @param genes A character vector with ensembl gene ids as 
#' 
#' @return Dataframe consisting of three character vectors of ensembl gene id, symbol and chromosome



gene_name_converter <- function(genes)
  
{
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

  ensembl_ids <- sub("[.][0-9]*","", genes)
  gene_ids <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol","chromosome_name"), values = ensembl_ids, mart = mart)
 
  
  return(gene_ids)
  
}

variance_stabilization <- function(edat, pdat, n_sv){
  sv_char <- rep("sv",n_sv)
  nsv_char <- as.character(seq(1,n_sv,1))
  sv_char_comp <- paste0(sv_char,nsv_char)
  plus_char <- c(rep("+", n_sv -1),"")
  sv_char_plus <- paste0(sv_char_comp, plus_char)
  design_char <- paste0("~ case +  sex + age + age*sex + case*sex +",paste(sv_char_plus, collapse = " "), "+ sex *(",paste(sv_char_plus, collapse = " "), ")") 
  design_form <- as.formula(design_char)
  de_seq_obj <- DESeqDataSetFromMatrix(countData = edat, colData = pdat, design = design_form)
  de_seq_sf  <- estimateSizeFactors(de_seq_obj)
  normalized_counts <- counts(de_seq_sf,normalized = TRUE)
  vst_deseq_obj <- vst(de_seq_sf, blind = FALSE)
  vsd_mat   <- assay(vst_deseq_obj)
  vsd_ret <- as.data.frame(cbind(t(vsd_mat), pdat$case, pdat$sex, pdat$age))
  colnames(vsd_ret) <- c(rownames(vsd_mat), "case","sex","age")
    
  return(vsd_ret)
  }