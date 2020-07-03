##%######################################################%##
#                                                          #
####               cancer cell drug datil               ####
#                                                          #
##%######################################################%##

drug_info <- fread("./intogen_cancer_drivers_actionability/Drug_details.tsv", skip = 9, key = "Tumor_type")
tumor_drug_info <- drug_info[!"N/A"]

drug_gene_target_info <- fread("./intogen_cancer_drivers_actionability/Protein_Drug_Interactions.tsv", skip = 9)


##%######################################################%##
#                                                          #
####                TCGA breast cancer                  ####
#                                                          #
##%######################################################%##

#load("./TCGA_RNA_seq.RData") #tumor_data
#load("./TCGA_RNA_Seq_normal_data.RData") #normal_data
#TCGA_sample_info <- fread("./sample_TCGA_RNA_Seq.csv")

BRCA_expr <- fread("./TCGA_BRCA_2014/gdac.broadinstitute.org_BRCA.mRNAseq_Preprocess.Level_3.2014101700.0.0/BRCA.uncv2.mRNAseq_RSEM_normalized_log2.txt")
selected_patient_ID <- unlist(lapply(colnames(BRCA_expr)[grep("11", colnames(BRCA_expr))], function(x) str_c(str_split(x, "-")[[1]][1:3], collapse = "-")))
#selected_sample <- colnames(BRCA_expr)[unlist(lapply(selected_patient_ID, function(x){y <- grep(x,colnames(BRCA_expr)); ifelse(length(y) >= 2, y, NA)}))]
selected_sample <- colnames(BRCA_expr)[unlist(lapply(selected_patient_ID, function(x){y <- grep(x,colnames(BRCA_expr)); if(length(y) >= 2) y}))]

filter_number <- unlist(lapply(selected_sample, function(x){y <- str_split(x, "-")[[1]][4]; ifelse(y=="11" | y=="01", "1", "2")}))
selected_sample <- selected_sample[-grep("2",filter_number)]
BRCA_expr <- BRCA_expr[, c("gene", selected_sample), with = F]
BRCA_expr <- as.data.frame(BRCA_expr)
rownames.gene <- unlist(lapply(BRCA_expr[, "gene"], function(x)str_split(x, "\\|")[[1]][2]))
#rownames.filter <- unlist(lapply(rownames.gene, function(x)ifelse(x=="?", "1", "2")))
#ignore_gene <- grep("1", rownames.filter)
#BRCA_expr <- BRCA_expr[-ignore_gene,]
rownames(BRCA_expr) <- rownames.gene
BRCA_expr <- BRCA_expr[,-1]

paired_info <- unlist(lapply(colnames(BRCA_expr), function(x)str_c(str_split(x, "-")[[1]][1:3], collapse = "-")))
group_info <- unlist(lapply(colnames(BRCA_expr), function(x){y <- str_split(x, "-")[[1]][4]; ifelse(y == "01", "T", "N")}))

source("./Tumor_Normal_function.R")
library("limma")

BRCA.foldchange <- get_paired_foldchange(BRCA_expr, group_info, paired_info)


##%######################################################%##
#                                                          #
####                        GSEA                        ####
#                                                          #
##%######################################################%##
source("./liver_cancer_function.R")

expres_gct(BRCA_expr, "BRCA")
phenotype_cls2(group_info, "BRCA")
