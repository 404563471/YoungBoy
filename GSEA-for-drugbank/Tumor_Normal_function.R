#group_info: "T" "T" "T" "N" "N" "N" "N"
#pair_info:"TCGA-E2-A1LS" "TCGA-BH-A0B3" "TCGA-E2-A1LH" "TCGA-BH-A1EW"
#length(group_info)==length(pair_info)

get_paired_foldchange <- function(expression, group_info, pair_info, pval = 0.01){
  expression <- na.omit(expression)
  paired <- factor(pair_info)
  tumor.v.normal <- factor(group_info, levels = c("N", "T"))# T - N
  design <- model.matrix(~paired+tumor.v.normal)
  fit <- lmFit(expression, design)
  fit <- eBayes(fit)
  fold_change <- topTable(fit, coef = "tumor.v.normalT", number = nrow(fit))
  fold_change <- fold_change[fold_change$adj.P.Val < pval,]
  return(fold_change)
}




transform_probeID <- function(gene_symbol, name){
              require(annotate)
              require(hgu133a.db)
              probeID <- select(hgu133a.db, keys = gene_symbol, keytype = "SYMBOL", columns = "PROBEID")
              probeID <- na.omit(probeID$PROBEID)
              probeID <- unique(probeID)
              write(probeID, file = paste0("~/cmap_TNBC_kinase/tag_list/", name, ".grp"))
              return(probeID)
}

str_right2.0 <- function(x, start_nu, stop_nu){
  substr(x, start_nu+nchar(x)+1, stop_nu+nchar(x)+1)
}


get_hc <- function(expression, method="com"){
  expression_cor <- cor(expression)
  expression_hc <- hclust(as.dist(1-expression_cor), method)
  return(expression_hc)
}


get_color_hc <- function(expression, name, group, estimate_score = NULL, method="com"){
  require("dendextend")
  color <- data.frame(color = rep(NA, length(colnames(expression))), stringsAsFactors = F)
  rownames(color) <- colnames(expression)
  
  if(!is.null(group)) {group_level <- levels(as.factor(group))}
  if(!is.null(estimate_score)){
    min_sturation <- min(estimate_score["TumorPurity",]) - 0.15
    multiple <- 1/(max(estimate_score["TumorPurity",])-min_sturation)
  }
  
  for(i in 1:length(group_level)){
    if(!is.null(estimate_score)){
      for(j in grep(group_level[i], group)){
        color[j, "color"] <- rainbow(length(group_level), s = (estimate_score["TumorPurity", j]-min_sturation)*multiple)[i]
      }
    }  
    else {
      color[grep(group_level[i], group), "color"] <- rainbow(length(group_level))[i]
    }
  }
  
  qc_hc <- get_hc(expression, method)
  hc_forplot <- as.dendrogram(qc_hc)
  labels_colors(hc_forplot) <- color$color[order.dendrogram(hc_forplot)]
  par(mai = c(2,0.1,0.1,0.1), fin = c(10, 6.99999))
  
  jpeg(paste("~/cmap_TNBC_kinase/hc_plot/", name, ".jpeg", sep = ""), width = 1500, height = 1000)
  plot(hc_forplot)
  dev.off()
}


##%######################################################%##
#                                                          #
####             estimate for tumor purity              ####
#                                                          #
##%######################################################%##
#how to download:
#install.packages("estimate", repos="http://R-Forge.R-project.org")

make_estimate <- function(expression, name){
  require("estimate")
  Document <- paste("~/cmap_TNBC_kinase/estimate/", name, "/", sep = "")
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