#sex_prediction algorithms with msbb datasets
# start with a dataframe consisting of sex, age and case 
# columns binded to a sample X gene matrix

source("~/RNA_SEQ/rna_seq_functions/MSBB_preprocess.R")
source("~/rna_seq_functions/rna_seq_functions.R")
library(e1071)
library(cluster)

gene_chr_table <- gene_name_converter(colnames(sva_df))
xy_chr_table <- gene_chr_table[gene_chr_table$chromosome_name == "Y" | gene_chr_table$chromosome_name == "X",]

xist_kd5md <- xy_chr_table[xy_chr_table$hgnc_symbol
                           %in% c("XIST","RPS4Y1"),]
colnames(sva_df) <-  sub("[.][0-9]*","", colnames(sva_df))
xy_df <- sva_df[,colnames(sva_df) %in% c(xy_chr_table$ensembl_gene_id,"sex")]
xist_df <- sva_df[,colnames(sva_df) %in% c(xist_kd5md$ensembl_gene_id,"sex")]
xy_df$sex <- as.factor(xy_df$sex)
set.seed(1)

#use 70% of dataset as training set and 30% as test set
sample <- sample(c(TRUE, FALSE), nrow(xy_df), replace=TRUE, prob=c(0.7,0.3))
train  <- xy_df[sample, ]
test   <- xy_df[!sample, ]


log_reg_model <- glm(sex ~ ., data = xy_df, family = binomial())
 tune.out <- tune(svm, sex~., data = train, kernel = "linear",
                             ranges = list(cost = c(0.1,1,10,100,1000)
                                           ))
 
 tune.best <- tune.out$best.model
 
 
 sex_prediction <- predict(tune.best, test)
 
 cross_tab_pred <- table(sex_prediction, test$sex)
 
 full_test <- predict(tune.best, xy_df)
 cross_tab_full <- table(full_test, xy_df$sex)
 
n_cols <- ncol(xy_df) - 1
 pam_full <- pam(xy_df[,1:n_cols],2, metric = "euclidean")
 
xy_df_pam <- as.data.frame(cbind(xy_df, pam_full$clustering))