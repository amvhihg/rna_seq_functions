pheno_data_t$case <- as.factor(pheno_data_t$case)
pheno_data_t$case <- relevel(pheno_data_t$case, ref = "Normal")
pheno_data_t <- pheno_data_t[order(rownames(pheno_data_t)),]
mayo_gene_cer <- mayo_gene_cer[,order(colnames(mayo_gene_cer))]
deseq_obj_sva <- DESeqDataSetFromMatrix(mayo_gene_cer, pheno_data_t, ~ case + sex + case*sex + age)
deseq_obj_sva <- estimateSizeFactors(deseq_obj_sva)
norm_cts <- counts(deseq_obj_sva, normalized = TRUE)

sva_df <- as.data.frame(cbind(t(norm_cts), pheno_data_t[,c("case","sex","age")]))

sva_obj <- sva_cleaning_mod(sva_df)

sva_obj$sv <- apply(sva_obj$sv,2,function(x){ scale(x, center = TRUE)})
sva_ret <- sva_obj$sv
sv_char <- rep("sv",sva_obj$n.sv)
nsv_char <- as.character(seq(1,sva_obj$n.sv,1))
sv_char_comp <- paste0(sv_char,nsv_char)

colnames(sva_ret) <- sv_char_comp
pheno_data_t <- as.data.frame(cbind(pheno_data_t, sva_ret))
pheno_data_t$age <- scale(pheno_data_t$age)
pheno_data_t$age <- as.double(pheno_data_t$age)


#DE_seq_obj        <- DESeqDataSetFromMatrix(countData = rosmap_final_count, colData = pheno_data_t, design = ~ case + sex + case:sex + sv1 +sv2+sv3 +(sv1+sv2+sv3)  *sex+ age + age:sex )
#DE_seq_obj        <- estimateSizeFactors(DE_seq_obj)
#normalized_counts <- counts(DE_seq_obj, normalized = TRUE)
#vst_de_seq_obj    <- vst(DE_seq_obj, blind = FALSE)
#vsd_mat           <- assay(vst_de_seq_obj)
#vsd_df            <- as.data.frame(cbind(t(vsd_mat), pheno_data_t$case,  pheno_data_t$sex,pheno_data_t$age))
#colnames(vsd_df) <- c(rownames(vsd_mat), "case","sex","age")
vsd_df <- variance_stabilization(edat = mayo_gene_cer, pdat = pheno_data_t, n_sv =sva_obj$n.sv)
pca_vsd_case <- pca_data_out(vsd_df, "case")
pca_vsd_sex  <- pca_data_out(vsd_df, "sex")
#pca_vsd_race <- pca_data_out(vsd_df, "race")
outliers_df <- pca_vsd_case[[2]]
outliers_df <- outliers_df$x
x_limit <- pca_vsd_case[[3]][[1]]
x_mean  <- pca_vsd_case[[3]][[3]]
x_point <- x_mean  - 3 * x_limit
x_point_r <- x_mean + 3 * x_limit
y_limit <- pca_vsd_case[[3]][[2]]
y_mean <- pca_vsd_case[[3]][[4]]
y_point <- y_mean - 3*y_limit
y_point_r <- y_mean + 3*y_limit
outliers_x_l <- names(which(outliers_df[,1] < x_point))
outliers_x_r <- names(which(outliers_df[,1] > x_point_r))
outliers_y <- names(which(outliers_df[,2] < y_point))
outliers_y_r <- names(which(outliers_df[,2] > y_point_r))
outliers <- c(outliers_x_l, outliers_y,outliers_x_r,outliers_y_r)
pheno_data_t <- pheno_data_t[!(rownames(pheno_data_t) %in% outliers),]
rosmap_df <- mayo_gene_cer[,!(colnames(mayo_gene_cer) %in% outliers)]
rownames(pheno_data_t) <- pheno_data_t$specimenID
rosmap_df <- rosmap_df[,colnames(rosmap_df) %in% rownames(pheno_data_t)]

#rosmap_df <- rosmap_final_count #for msbb datasets
pheno_data_t$case <- as.factor(pheno_data_t$case)
pheno_data_t$case <- relevel(pheno_data_t$case, ref = "Normal")
pheno_data_t$int <- ifelse(pheno_data_t$sex == "male" & pheno_data_t$case == "AD","B","A")
pheno_data_t$int <- as.factor(pheno_data_t$int)
rownames(rosmap_df) <-   sub("[.][0-9]*","", rownames(rosmap_df))
full_model_list <- de_seq_function(rosmap_df, pheno_data_t, sva_obj$n.sv)
full_model   <- full_model_list[[1]]
size_factors <- full_model_list[[2]]



fem_deqseq <- size_factors[,size_factors$sex =="female"]
male_deseq <- size_factors[,size_factors$sex == "male"]
fem_model <- sex_strat_deseq(fem_deqseq, sva_obj$n.sv)

male_model <- sex_strat_deseq(male_deseq, sva_obj$n.sv)
full_model <- DESeq(size_factors)


results_int <- results(full_model,contrast = c("int","B","A"),independentFiltering = FALSE)
int_table <- do.call("cbind",results_int@listData)
int_table <- as.data.frame(int_table)
rownames(int_table) <- results_int@rownames

results_case <- results(full_model,contrast = c("case","AD", "Normal"),independentFiltering = FALSE)
case_table   <- do.call("cbind",results_case@listData)
case_table <- as.data.frame(case_table)
rownames(case_table) <- results_case@rownames

results_sex <- results(full_model, contrast = c("sex","male","female"), independentFiltering = FALSE)
sex_table <- do.call("cbind",results_sex@listData)
sex_table <- as.data.frame(sex_table)
rownames(sex_table) <- results_sex@rownames


case_table_sym <- gene_name_converter(rownames(case_table))
case_table_sym$dup <- duplicated(case_table_sym$ensembl_gene_id)
case_table_sym <- case_table_sym[case_table_sym$dup == FALSE,]
rownames(case_table_sym) <- case_table_sym$ensembl_gene_id
case_table_fin <- merge(case_table, case_table_sym, by = "row.names")

case_table_fin <- case_table_fin[,2:ncol(case_table_fin)]

sex_table_fin <- merge(sex_table, case_table_sym, by = "row.names")
int_table_fin <- merge(int_table, case_table_sym, by = "row.names")
fem_model <- DESeq(fem_deqseq)
male_model <- DESeq(male_deseq)
results_fem <- results(fem_model,contrast = c("case","AD" ,"Normal"),independentFiltering =FALSE)
fem_table   <- do.call("cbind",results_fem@listData) 
fem_table   <- as.data.frame(fem_table)
rownames(case_table) <- results_case@rownames
rownames(fem_table) <- results_fem@rownames

results_male <- results(male_model, contrast = c("case", "AD","Normal"), independentFiltering = FALSE)
male_table <- do.call("cbind", results_male@listData)
male_table <- as.data.frame(male_table)
rownames(male_table) <- results_male@rownames
#rownames(male_table) <-  sub("[.][0-9]*","", rownames(male_table))


fem_table_sym <- gene_name_converter(rownames(fem_table))

fem_table_sym$dup <- duplicated(fem_table_sym$ensembl_gene_id)
fem_table_sym <- fem_table_sym[fem_table_sym$dup   == FALSE,]


rownames(fem_table_sym) <- fem_table_sym$ensembl_gene_id



fem_table_fin <- merge(fem_table, fem_table_sym, by  = "row.names")

rownames(fem_table_fin) <- fem_table_fin$Row.names
fem_table_fin <- fem_table_fin[,2:ncol(fem_table_fin)]
case_fem_comp <- merge( case_table_fin, fem_table_fin, by = "row.names")
head(case_table)
male_table_fin <- merge(male_table, fem_table_sym, by = "row.names")

int_table_fin <- merge(int_table, fem_table_sym, by = "row.names")



sex_table_fin <- merge(sex_table, fem_table_sym, by = "row.names")

list_tables <- list(sex_table_fin, int_table_fin, case_table_fin, fem_table_fin, male_table_fin)

msbb_bm22_nofilt_master <- do.call("cbind",list_tables)
save.image("MSBB_noindfilt_version.RData")
write.csv(msbb_bm22_nofilt_master, "msbb_bm22_master_sign_cor.csv")

write.csv(int_table_fin, "int_table_nofilt.csv")
write.csv(male_table_fin, "male_nofilt.csv")
write.csv(fem_table_fin, "fem_nofilt.csv")