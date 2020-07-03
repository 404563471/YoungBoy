source("~/liver_cancer/liver_cancer_function.R")

library(affy)
library(foreach)
library(doParallel)

setwd("~/liver_cancer/test/")
test_raw <- ReadAffy()


cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
qc_result <- foreach(GSE_names = list("GSE40873", "GSE29721", "GSE19665", "GSE36076", "GSE58208", "GSE49515", "GSE62232", "GSE6764", "GSE17548", "GSE9843")) %dopar% {
              file_name <- paste("~/liver_cancer/", GSE_names, "/", sep = "")
              setwd(file_name)
              rawdata <- ReadAffy()
              make_qc2.0(rawdata)
}
stopCluster(cl)

GSE_names <- c("GSE40873", "GSE29721", "GSE19665", "GSE36076", "GSE58208", "GSE49515", "GSE62232", "GSE6764", "GSE17548", "GSE9843")
for(i in 1:length(GSE_names)){
  
  name <- paste(GSE_names[i], "qc_result", sep = "_")
  assign(name, qc_result[[i]])
  print(name)
  print(get(name))
  
}



setwd("~/liver_cancer/GSE9843/")
GSE9843_raw <- ReadAffy()

GSE9843_qc_result <- make_qc2.0(GSE9843_raw)



library(bigmemory)
#create a matrix of size 1GB aproximatelly
n <- 10000
m <- 10000
c <- matrix(runif(n*m),n,m)
#convert it to bigmatrix
x <- as.big.matrix(x = c, type = "double", 
                   separated = FALSE, 
                   backingfile = "example.bin", 
                   descriptorfile = "example.desc")
# get a description of the matrix
mdesc <- describe(x)
# Create the required connections    
cl <- makeCluster(detectCores ())
registerDoSNOW(cl)
## 1) No referencing
out <- foreach(linID = 1:4, .combine=c) %dopar% {
  t <- attach.big.matrix("example.desc")
  for (i in seq_len(30L)) {
    for (j in seq_len(m)) {
      y <- t[i,j]
    }
  }
  return(0L)
}


e1 <- new.env()
library(Rdsm)



###########################################################
########## after QC
##########
###########################################################
probe_gene <- data.frame(list(probe = c("229090_at", "224945_at", "215121_x_at", "238619_at", "230930_at", "1563298_at", "239503_at"), gene = c("ZEB1_AS1", "ENST00000417262", "NR_033661", "AK129699", "BC043009", "ENST00000449772", "BC041955")), stringsAsFactors = F)
probe_gene[8, ] <- c("232230_at", "OLMALINC")
probe_gene[9:10,] <- c("200831_s_at", "200832_s_at", "SCD1", "SCD2")
probe_gene[11:17,] <- c("210795_s_at", "217797_at", "217494_s_at", "239153_at", "212495_at", "210662_at", "1561542_at", "MEG3", "UFC1", "PTENP1", "HOTAIR", "BC015378", "UC002TVK", "AK056988")
probe_gene[18:19,] <- c("1568611_at", "225457_s_at", "BC013423", "OLMALINC")



GSE40873_cel <- list.celfiles("~/liver_cancer/GSE40873/")
GSE40873_cel <- intersect(GSE40873_cel, GSE40873_qc_result$good_sample)
GSE40873_expression <- get_expression1.0("GSE40873", GSE40873_cel)
GSE40873_clinical <- read.delim("~/liver_cancer/GSE40873/GSE40873_series_matrix.txt", header = T, stringsAsFactors = F, sep = ",")#11 rows is survial 
colnames(GSE40873_clinical) <- GSE40873_clinical[1,]
#GSE40873_expression <- lapply(as.list(GSE40873_clinical[1,]), FUN = function(x){y <- grep(x, colnames(GSE40873_expression)); GSE40873_expression["survival_day", y] <- GSE40873_clinical[11,x]; GSE40873_expression["events", y] <- GSE40873_clinical[12,x]})
GSE40873_clinical_expression <- GSE40873_expression


GSE40873_expression_gcrma <- get_expression1.0("GSE40873", "gcrma", not_qc = F)
GSE40873_clinical_expression <- GSE40873_expression_gcrma

GSE40873_expression_rma <- get_expression1.0("GSE40873", "rma", GSE40873_cel)
GSE40873_clinical_expression <- GSE40873_expression_rma


GSE40873_expression_mas5 <- get_expression1.0("GSE40873", "mas5", GSE40873_cel)
GSE40873_clinical_expression <- GSE40873_expression_mas5

for(i in GSE40873_clinical[1,]){
  y <- grep(i, colnames(GSE40873_clinical_expression))
  if (length(y) == 0) next()
  else {GSE40873_clinical_expression["survival_day", y] <- GSE40873_clinical[11,i]
        GSE40873_clinical_expression["events", y] <- GSE40873_clinical[12,i]}
  # i don't know why after this operation, nueerics of data.fream become character........ maybe everycol have to be one class
}
GSE40873_clinical_expression["survival_day",] <- lapply(GSE40873_clinical_expression["survival_day",], FUN = function(x){strsplit(x, ": ")[[1]][2]})
GSE40873_clinical_expression["events",] <- lapply(GSE40873_clinical_expression["events",], FUN = function(x){strsplit(x, ": ")[[1]][2]})
GSE40873_clinical_expression["MO_status",] <- grepl("^multicentric occurrence", GSE40873_clinical_expression["events",])

LincRNA_survival <- as.data.frame(t(GSE40873_clinical_expression[c("survival_day", probe_gene$probe),]))
LincRNA_survival <- as.data.frame(apply(LincRNA_survival, 2, FUN = function(x){as.numeric(x)}))
LincRNA_survival$MO_status <- t(GSE40873_clinical_expression["MO_status",])
LincRNA_survival$MO_status <- as.logical(t(GSE40873_clinical_expression["MO_status",]))

library(survminer)
library(survival)
library(progress)
pb <- progress_bar$new(format = "survival ploting: (:spin) [:bar] :percent", total = nrow(probe_gene), clear = FALSE, width = 100)
for (i in 1:nrow(probe_gene)) {
      y <- Surv(LincRNA_survival$survival_day, LincRNA_survival$MO_status)
      probe_cutoff <- cutoff.survival(LincRNA_survival, probe_gene$probe[i])
      if (probe_cutoff$pvalue == 1) print(probe_gene$gene[i])
      else 
      group <- LincRNA_survival[, probe_gene$probe[i]] <= probe_cutoff$cutoff
      kmfit <- survfit(y~group)
      jpeg(paste("~/liver_cancer/GSE40873_survival/GSE40873", probe_gene$gene[i], Sys.Date(), "mas5", ".jpeg", sep = "-")) # 2017-01-03
      print(ggsurvplot(kmfit, pval = T, legend.title = paste("GSE40873-", probe_gene$gene[i], sep = ""), legend.labs = c("high_expression", "low_expression"), risk.table = T))
      dev.off()
      pb$tick()
}

method <- c("rma", "gcrma", "mas5")
GSE40873_old_expression <- get_expression2.0("GSE40873")
  for(k in 1:length(method)){
    GSE_expression <- paste("GSE40873_old_expression", method[k], sep = "_")
    assign(GSE_expression, GSE40873_old_expression[[k]])
    GSE40873_clinical_expression <- get(GSE_expression)
    
    for(i in GSE40873_clinical[1,]){
      y <- grep(i, colnames(GSE40873_clinical_expression))
      if (length(y) == 0) next()
      else {GSE40873_clinical_expression["survival_day", y] <- GSE40873_clinical[11,i]
            GSE40873_clinical_expression["events", y] <- GSE40873_clinical[12,i]}
      # i don't know why after this operation, nueerics of data.fream become character........ maybe everycol have to be one class
    }
    GSE40873_clinical_expression["survival_day",] <- lapply(GSE40873_clinical_expression["survival_day",], FUN = function(x){strsplit(x, ": ")[[1]][2]})
    GSE40873_clinical_expression["events",] <- lapply(GSE40873_clinical_expression["events",], FUN = function(x){strsplit(x, ": ")[[1]][2]})
    GSE40873_clinical_expression["MO_status",] <- grepl("^multicentric occurrence", GSE40873_clinical_expression["events",])
    
    LincRNA_survival <- as.data.frame(t(GSE40873_clinical_expression[c("survival_day", probe_gene$probe),]))
    LincRNA_survival <- as.data.frame(apply(LincRNA_survival, 2, FUN = function(x){as.numeric(x)}))
    LincRNA_survival$MO_status <- t(GSE40873_clinical_expression["MO_status",])
    LincRNA_survival$MO_status <- as.logical(t(GSE40873_clinical_expression["MO_status",]))
    for (i in 1:nrow(probe_gene)) {
          y <- Surv(LincRNA_survival$survival_day, LincRNA_survival$MO_status)
          probe_cutoff <- cutoff.survival(LincRNA_survival, probe_gene$probe[i])
          if (probe_cutoff$pvalue == 1) print(probe_gene$gene[i])
          else 
          group <- LincRNA_survival[, probe_gene$probe[i]] <= probe_cutoff$cutoff
          kmfit <- survfit(y~group)
          jpeg(paste("~/liver_cancer/GSE40873_survival/GSE40873_old", probe_gene$gene[i], Sys.Date(), method[k], ".jpeg", sep = "-")) # 2017-01-03
          print(ggsurvplot(kmfit, pval = T, legend.title = paste("GSE40873-", probe_gene$gene[i], sep = ""), legend.labs = c("high_expression", "low_expression"), risk.table = T))
          dev.off()
          print(i)
    }
  }

#save(LincRNA_survival, file = "~/liver_cancer/LincRNA_survival.RData")

y <- Surv(LincRNA_survival$survival_day, LincRNA_survival$MO_status)
probe_cutoff <- cutoff.survival(LincRNA_survival, "1568611_at")
group <- LincRNA_survival[, "1568611_at"] <= probe_cutoff$cutoff
kmfit <- survfit(y~group)
jpeg(paste("~/liver_cancer/BC013423-survival", Sys.Date(), ".jpeg", sep = "-"), quality = 100, width = 600, height = 800)
print(ggsurvplot(kmfit, size = 0.8, pval = T, surv.scale = "percent", ylim = c(0.2, 1), legend.title = NULL, legend = c(0.8,0.9),
                 legend.labs = c("lncRNA-PE(+)","lncRNA-PE(-)"), risk.table = T, xlab = "Postoperative days", 
                 ylab = "Multicentric Occurrence Free Survival rate", risk.table.col = "strata", risk.table.title = NULL, 
                 risk.table.y.text.col = F, risk.table.fontsize = 6, font.main = "red", font.x = c(18, "plain", "black"),
                 font.y = c(18, "plain", "black"), font.tickslab = c(14, "plain", "black"), font.legend = c(16, "plain", "black")))
dev.off()


#################################################
####GSE40873_GSEA
#################################################
GSE40873_unique_expression_rma <- get_unique_expression(GSE40873_expression_rma, "hgu133plus2.db")
GSE40873_entrezg_expression_rma <- get_expression1.0("GSE40873", "rma", custom_CDF = "entrezg", not_qc = F)

make_all_GSEA_files(list(GSE40873_expression_rma, GSE40873_unique_expression_rma, GSE40873_entrezg_expression_rma), "GSE40873", probe_gene)



##############################################################################
## GSE29721 ###
##############################################################################
GSE29721_matrix <- read.delim("~/liver_cancer/GSE29721/GSE29721_series_matrix.txt", stringsAsFactors = F, row.names = 1)

GSE29721_cel <- list.celfiles("~/liver_cancer/GSE29721/")
GSE29721_cel <- intersect(GSE29721_cel, GSE29721_qc_result$good_sample)
GSE29721_expression_rma <- get_expression1.0("GSE29721", "rma", GSE29721_cel)

GSE29721_qc_number <- lapply(colnames(GSE29721_matrix), FUN = function(x){grep(x, colnames(GSE29721_expression_rma))})
GSE29721_qc_group <- c()
for(i in 1:length(GSE29721_qc_number)){
  number <- GSE29721_qc_number[[i]]
  if(length(number) == 0) next()
  else GSE29721_qc_group[number] <- GSE29721_matrix[1,i]
}

GSE29721_qc_group <- sub("HCC tissue", "Tumor", GSE29721_qc_group)
GSE29721_qc_group <- sub("Normal liver tissue", "Normal", GSE29721_qc_group)

GSE29721_dif <- get_flodchange(GSE29721_qc_group, GSE29721_expression_rma)
for(i in 1:nrow(probe_gene)){
 jpeg(paste("~/liver_cancer/dif_plot/GSE29721/", probe_gene$gene[i], "_dif_", Sys.Date(),".jpeg", sep = ""))     
 print(get_plot(probe_gene$probe[i], probe_gene$gene[i], GSE29721_expression_rma, GSE29721_qc_group, GSE29721_dif, "GSE29721"))
 dev.off()
 print(i)
} 

GSE29721_origin_expression_rma <- get_expression1.0("GSE29721", "rma")
GSE29721_group <- get_group(GSE29721_matrix, GSE29721_origin_expression_rma, 1, c("HCC tissue", "Normal liver tissue"))


GSE29721_origin_dif <- get_flodchange(GSE29721_group, GSE29721_origin_expression_rma)
for(i in 1:nrow(probe_gene)){
 jpeg(paste("~/liver_cancer/dif_plot/GSE29721/", probe_gene$gene[i], "_origin_dif_", Sys.Date(),".jpeg", sep = ""))     
 print(get_plot(probe_gene$probe[i], probe_gene$gene[i], GSE29721_origin_expression_rma, GSE29721_group, GSE29721_origin_dif, "GSE29721"))
 dev.off()
 print(i)
} 


GSE29721_origin_unique_expression <- get_unique_expression(GSE29721_origin_expression_rma, "hgu133plus2.db")
GSE29721_unique_expression <- get_unique_expression(GSE29721_expression_rma, "hgu133plus2.db")

GSE29721_origin_entrezg_expression <- get_expression1.0(GSE_name = "GSE29721", method = "rma", custom_CDF = "entrezg")
GSE29721_entrezg_expression <- get_expression1.0("GSE29721", "rma", "entrezg", not_qc = F)

make_all_GSEA_files(list(GSE29721_expression_rma, GSE29721_unique_expression, GSE29721_entrezg_expression), "GSE29721", probe_gene, sample_group = GSE29721_qc_group)
##############################################################
#### {b} Blood GSE36076 GSE58028 GSE49515 #####################
##############################################################

cl <- makeCluster(3, type = "FORK")
registerDoParallel(cl)
blood_rma_expression <- foreach(i = list("GSE36076", "GSE58208", "GSE49515")) %dopar% {
        get_expression1.0(GSE_name = i, method = "rma", not_qc = F)
}
stopImplicitCluster()

blood_name <- c("GSE36076", "GSE58208", "GSE49515")
for (i in 1:length(blood_name)){
GSE_expression <- paste(blood_name[i], "_expression", sep = "")
assign(GSE_expression, blood_rma_expression[[i]])
delete_sample <- Reduce(union, list(grep("PAN", colnames(get(GSE_expression))), grep("GAS", colnames(get(GSE_expression))), grep("CHB", colnames(get(GSE_expression)))))
assign(GSE_expression, get(GSE_expression)[, -delete_sample])

GSE_group <- paste(blood_name[i], "group", sep = "_")
assign(GSE_group, grepl("HCC", colnames(get(GSE_expression))))
assign(GSE_group, sub(TRUE, "Tumor", get(GSE_group)))
assign(GSE_group, sub(FALSE, "Normal", get(GSE_group)))
}


for(i in blood_name){
  GSE_expression <- paste(i, "expression", sep = "_")
  GSE_group <- paste(i, "group", sep = "_")
  GSE_dif <- paste(i, "dif", sep = "_")
  
  assign(GSE_dif, get_flodchange(get(GSE_group), get(GSE_expression)))
  for(j in 1:nrow(probe_gene)){
    GSE_plot <- paste(i, probe_gene$gene, sep = "_")
    assign(GSE_plot, get_plot(probe_gene$probe[j], probe_gene$gene[j], get(GSE_expression), get(GSE_group), get(GSE_dif), i))
    ggsave(paste("~/liver_cancer/dif_plot/blood/", i, "_", probe_gene$gene[j], "_", Sys.Date(), ".jpeg", sep = ""), plot = get(GSE_plot), device = "jpeg")
    cat(i, j, sep = "\n")
  }
  
}

for(i in blood_name){
  GSE_expression <- paste(i, "expression", sep = "_")
  GSE_hc <- paste(i, "hc", sep = "_")
  assign(GSE_hc, get_qc_hc(get(GSE_expression)))
  print(plot(get(GSE_hc)))
}

library(sva)
blood_expression <- as.matrix(data.frame(GSE36076_expression, GSE58208_expression))
csif <- data.frame("sample" = colnames(blood_expression), "batch" = rep(c(1,2), c(ncol(GSE36076_expression), ncol(GSE58208_expression))), "condition" = c(GSE36076_group, GSE58208_group))
modcombat <- model.matrix(~1, data = csif)
batch <- csif$batch
blood_expression <- ComBat(dat = blood_expression, batch = batch, mod = modcombat, par.prior = T, prior.plots = T)
blood_hc <- get_qc_hc(blood_expression)
plot(blood_hc)

blood_expression <- as.data.frame(blood_expression)
blood_group <- c(GSE36076_group, GSE58208_group)
blood_dif <- get_flodchange(blood_group, blood_expression)
for(i in 1:nrow(probe_gene)){
  GSE_plot <- paste("Blood", probe_gene$gene[i], sep = "_")
  assign(GSE_plot, get_plot(probe_gene$probe[i], probe_gene$gene[i], blood_expression, blood_group, blood_dif, "Blood"))
  ggsave(paste("~/liver_cancer/dif_plot/blood/", GSE_plot, "_", Sys.Date(), ".jpeg", sep = ""), plot = get(GSE_plot))
}



GSE36076_old_expression <- get_expression1.0("GSE36076", "rma")
delete_sample <- Reduce(union, list(grep("PAN", colnames(GSE36076_old_expression)), grep("GAS", colnames(GSE36076_old_expression)), grep("CHB", colnames(GSE36076_old_expression))))
GSE36076_old_expression <- GSE36076_old_expression[, -delete_sample]
GSE36076_old_group <-  grepl("HCC", colnames(GSE36076_old_expression))
GSE36076_old_group <-  sub(TRUE, "Tumor", GSE36076_old_group)
GSE36076_old_group <- sub(FALSE, "Normal", GSE36076_old_group)

GSE36076_old_dif <- get_flodchange(GSE36076_old_group, GSE36076_old_expression)
for(i in 1:nrow(probe_gene)){
  GSE_plot <- paste("GSE36076_old", probe_gene$gene[i], sep = "_")
  assign(GSE_plot, get_plot(probe_gene$probe[i], probe_gene$gene[i], GSE36076_old_expression, GSE36076_old_group, GSE36076_old_dif, "GSE36076_old"))
  ggsave(paste("~/liver_cancer/dif_plot/blood/", GSE_plot, "_", Sys.Date(), ".jpeg", sep = ""), plot = get(GSE_plot), device = "jpeg")
}


GSE58208_origin_expression <- get_expression1.0("GSE58208", "rma")
GSE58208_origin_expression <- GSE58208_origin_expression[, -grep("HCC", colnames(GSE58208_origin_expression))]
GSE58208_origin_group <-  grepl("CHB", colnames(GSE58208_origin_expression))
GSE58208_origin_group <-  sub(TRUE, "Tumor", GSE58208_origin_group)
GSE58208_origin_group <- sub(FALSE, "Normal", GSE58208_origin_group)
GSE58208_origin_dif <- get_flodchange(GSE58208_origin_group, GSE58208_origin_expression)
for(i in 1:nrow(probe_gene)){
  GSE_plot <- paste("GSE58208_origin", probe_gene$gene[i], sep = "_")
  assign(GSE_plot, get_plot(probe_gene$probe[i], probe_gene$gene[i], GSE58208_origin_expression, GSE58208_origin_group, GSE58208_origin_dif, "GSE58208_origin"))
  ggsave(paste("~/liver_cancer/dif_plot/blood/", GSE_plot, "_", Sys.Date(), ".jpeg", sep = ""), plot = get(GSE_plot), device = "jpeg")
}


#################################################
##
##GSE84044
################################################
GSE84044_raw <- get_rawdata1.0("~/liver_cancer/GSE84044/")
GSE84044_qc_result <- make_qc2.0(GSE84044_raw)

GSE84044_matrix <- read.delim("~/liver_cancer/GSE84044/GSE84044_series_matrix.txt", header = T, stringsAsFactors = F)

GSE84044_expression_rma <- get_expression1.0("GSE84044", "rma", not_qc = F)
GSE84044_expression_rma_beforsva_hc <- get_qc_hc(GSE84044_expression_rma)
plot(GSE84044_expression_rma_beforsva_hc)

GSE84044_group <- get_group(GSE84044_matrix, GSE84044_expression_rma, 4)

GSE84044_expression_rma_sva <- get_combat_expression(GSE84044_expression_rma, GSE84044_group)
get_color_hc(GSE84044_expression_rma_sva, GSE84044_group, "GSE84044_after_sva")

make_estimate(GSE84044_expression_rma_sva, "GSE84044")

GSE84044_estimate <- read.delim("~/liver_cancer/estimate_result/GSE84044/GSE84044estimate_score.gct", row.names = 1)
GSE84044_estimate <- GSE84044_estimate[, -1]
get_color_hc2.0(GSE84044_expression_rma_sva, GSE84044_estimate, "GSE84044_estimate")

GSE84044_expression_rma_sva_unique <- get_unique_expression(GSE84044_expression_rma_sva)
GSE84044_expression_rma_sva_entrez <- get_expression1.0("GSE84044", "rma", "entrezg", not_qc = F)
make_all_GSEA_files(list(GSE84044_expression_rma_sva, GSE84044_expression_rma_sva_unique, GSE84044_expression_rma_sva_entrez), "GSE84044", probe_gene)




library(simpleaffy)
GSE84044_PMA <- detection.p.val(get_rawdata1.0("GSE84044", not_qc = F))
GSE84044_PMA_call <- GSE84044_PMA$call
#GSE84044_threshold <- apply(GSE84044_PMA_call[probe_gene$probe,], 1, FUN = function(x) sum(as.numeric(table(x, exclude = "A"))))
#GSE84044_threshold <- data.frame(GSE84044_threshold, probe_gene)
GSE84044_threshold <- get_present_threshold(GSE84044_PMA, probe_gene)
probe_gene <- probe_gene[-c(5,11,13,14,17),]
probe_gene <- probe_gene[-12,]
probe_gene$probe[2] <- "231233_at"


GSE84044_expression_rma_filter <- get_expression1.0("GSE84044", "rma", 43, GSE84044_PMA, not_qc = F)
GSE84044_expression_rma_filter_sva <- get_combat_expression(GSE84044_expression_rma_filter, batch = GSE84044_group)

GSE84044_cor <- get_cor2.0(GSE84044_expression_rma_filter_sva, "224945_at")
GSE84044_GO <- GO_analyse(GSE84044_expression_rma_filter_sva, GSE84044_cor, name = "GSE84044_rma_filter_entrez", pos_threshold = 0.7)
library(ggplot2)

GSE84044_GO_pos <- GSE84044_GO$pos_GO[GSE84044_GO$pos_GO$Size > 15 & GSE84044_GO$pos_GO$Size < 500, ]
GSE84044_GO_pos <- GSE84044_GO_pos[order(GSE84044_GO_pos$OddsRatio, decreasing = T),]

GSE84044_GO_plot <- ggplot(data = GSE84044_GO_pos, aes(x = Term, order = OddsRatio)) + 
  scale_y_continuous(trans = "log2") +  geom_bar(aes(y = Size), stat = "identity", show.legend = F) + 
  geom_bar(aes(y = Count, fill = Term, color = Term), stat = "identity", show.legend = F) +
  geom_bar(aes(y = ExpCount), stat = "identity", show.legend = F) +
  geom_text(aes(y = 300, label = Pvalue), position = "dodge", size = 10)+ coord_flip() + theme(axis.text = element_text(size = 30))

ggsave("~/liver_cancer/GOstats/test.jpeg", plot = GSE84044_GO_plot, device = "jpeg", width = 50, height = 30, limitsize = F)






GSE84044_expression_rma_filter_sva_entrez <- get_unique_expression(GSE84044_expression_rma_filter_sva, symbol = F)
GSE84044_expression_rma_filter_sva_symbol <- get_unique_expression(GSE84044_expression_rma_filter_sva, symbol = T)

make_all_GSEA_files(list(GSE84044_expression_rma_filter_sva, GSE84044_expression_rma_filter_sva_symbol, GSE84044_expression_rma_filter_sva_entrez),
                    "GSE84044", probe_gene)


####################################
###
###PVCA expample
####################################

library(golubEsets)
library(pvca)
data("Golub_Merge")
pct_threshold <- 0.6
batch.factors <- c("ALL.AML", "BM.PB", "Source")
pvcaObj <- pvcaBatchAssess(Golub_Merge, batch.factors, pct_threshold)
bp <- barplot(pvcaObj$dat, xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2, main="PVCA estimation bar chart")
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)



