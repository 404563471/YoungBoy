##%######################################################%##
#                                                          #
####  produce phenotype cls file for Tumor and Normal   ####
#                                                          #
##%######################################################%##
phenotype_cls2 <- function(group_info, name){
  filename <- paste0("~/GSEA-P-R/Datasets/", name, ".cls")
  first <- ifelse(group_info[1]=="T", "Tumor", "Normal")
  second <- ifelse(group_info[1]=="T", "Normal", "Tumor")
  write(paste(length(group_info), "2", "1", sep = "\t"),file = filename, sep = "\t")
  write(paste("#", first, second, sep = "\t"), file = filename, sep = "\t", append = T)
  write(group_info, file = filename, sep = "\t", append = T, ncolumns = length(group_info))
}

##%######################################################%##
#                                                          #
####            produce gene expression gct             ####
####       file, saved in "~/GSEA-P-R/Datasets/"        ####
#                                                          #
##%######################################################%##
expres_gct <- function(expression, name){
  expression <- na.omit(expression)
  expression_GSEA <- cbind(Description = "na", expression)
  filename_exprs <- paste("~/GSEA-P-R/Datasets/", name, ".gct", sep = "")
  write("#1.2", file = filename_exprs)
  write(as.character(c(nrow(expression), ncol(expression))), file = filename_exprs, ncolumns = 2, sep = "\t", append = T)
  write(c("NAME", colnames(expression_GSEA)), file = filename_exprs, ncolumns = ncol(expression_GSEA)+1, sep = "\t", append = T)
  write.table(expression_GSEA, file = filename_exprs, sep = "\t", append = T, quote = F, col.names = F, row.names = T)
}

##%######################################################%##
#                                                          #
####              excute gsea command line              ####
#                                                          #
##%######################################################%##

#java -cp gsea2-2.2.3.jar -Xmx2048m xtools.gsea.Gsea -res Datasets/BRCA.gct -cls Datasets/BRCA.cls -gmx GeneSetDatabases/Human_DrugBank_all_entrezgene.gmt -collapse false

##%######################################################%##
#                                                          #
####             estimate for tumor purity              ####
#                                                          #
##%######################################################%##
#how to download:
#install.packages("estimate", repos="http://R-Forge.R-project.org")

make_estimate <- function(expression, name){
  require("estimate")
  Document <- paste("~/R_work/GSEA-for-drugbank/", name, "/", sep = "")
  filename <- paste(name, "estimate", sep = "")
  dir.create(Document)
  write.table(expression, paste(Document, filename, ".txt", sep = ""), sep = "\t", quote = F)
  filterCommonGenes(input.f = paste(Document, filename, ".txt", sep = ""), output.f = paste(Document, filename, ".gct", sep = ""), id = "GeneSymbol")
  estimate_score_gct <- paste(Document, filename, "_score", ".gct", sep = "") 
  estimateScore(paste(Document, filename, ".gct", sep = ""), estimate_score_gct, platform="illumina")
  
  estimate_score <- read.delim(estimate_score_gct, row.names = 1, stringsAsFactors = F, skip = 1)
  estimate_score <- estimate_score[, -1]
  return(estimate_score)
}