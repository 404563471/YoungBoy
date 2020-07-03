library(affy) 
library(tcltk)
library(CLL)
library(annotate)
library(hgu133plus2.db)

#setwd("/home/yhy/liver_cancer/GSE45267/")
#cel.files <- list.files(pattern = ".cel$", ignore.case = T, full.names = T)
#GSE45267_rawdata <- read.affybatch(cel.files)
#GSE45267_gcrma <- gcrma::gcrma(GSE45267_rawdata)
#GSE45267_expression <- exprs(GSE45267_gcrma)
#GSE45267_expression <- as.data.frame(GSE45267_expression)

probe_symbol <- as.data.frame(getSYMBOL(rownames(GSE45267_expression), "hgu133plus2.db"))
colnames(probe_symbol) <- "SYMBOL"
GSE45267_expression <- GSE45267_expression[!is.na(getSYMBOL(rownames(GSE45267_expression), "hgu133plus2.db")), ]


GSE45267_group <- apply(as.array(colnames(GSE45267_expression)), 1, FUN = function(x)str_right(x,5))
GSE45267_expression[, c(grep("1", group), grep("2", group))] <- NULL
GSE45267_group <- apply(as.array(colnames(GSE45267_expression)), 1, FUN = function(x)str_right(x,5))
GSE45267_group <- sub("T", "Tumor", GSE45267_group)
GSE45267_group <- sub("N", "Normal", GSE45267_group)


GSE45267_dif <- get_flodchange(group, GSE45267_expression, 1)#1.29


get_cor <- function(eset_name, eset, group, T, met){eset_Tumor <- eset[, grep(T, group)] # "T" or "Tumor"
            cor <- apply(eset_Tumor, 1, FUN = function(x){cor(as.numeric(eset_Tumor["232230_at", ]), as.numeric(x), method = met,use = "pair")})
            # "pearson" (default), "kendall", or "spearman"
            cor <- as.data.frame(cor)
            GSE_cor <- paste(eset_name, "cor", sep = "_")
            colnames(cor) <- GSE_cor
            GSE_p <- paste(eset_name, "P.value", sep = "_")
            cor[, GSE_p] <- apply(eset_Tumor, 1, FUN = function(x){ y = cor.test(as.numeric(eset_Tumor["232230_at", ]), as.numeric(x), alternative = "two.sided", method = met);return(y$p.value)})
            cor$probe_ID <- rownames(cor)
            cor_0.3 <- cor[abs(cor[, GSE_cor]) > 0.3 & cor[, GSE_p] < 0.05, ]
            cor_0.3 <- na.omit(cor_0.3)
            cor_SYMBOL <- as.data.frame(getSYMBOL(rownames(cor_0.3), "hgu133plus2.db"))
            cor_SYMBOL$probe_ID <- rownames(cor_SYMBOL)
            cor_0.3 <- merge(cor_0.3, cor_SYMBOL, by.x = "probe_ID")
            colnames(cor_0.3)[4] <- "SYMBOL"
            #cor_0.3 <- cor_0.3[order(cor_0.3[, 2], decreasing = T), ]
            return(cor_0.3)
}
get_cor_2.0 <- function(eset_name, eset, group, T, met, probe){eset_Tumor <- eset[, grep(T, group)] # "T" or "Tumor"
            cor <- apply(eset_Tumor, 1, FUN = function(x){cor(as.numeric(eset_Tumor[probe, ]), as.numeric(x), method = met,use = "pair")})
            # "pearson" (default), "kendall", or "spearman"
            cor <- as.data.frame(cor)
            GSE_cor <- paste(eset_name, "cor", sep = "_")
            colnames(cor) <- GSE_cor
            GSE_p <- paste(eset_name, "P.value", sep = "_")
            cor[, GSE_p] <- apply(eset_Tumor, 1, FUN = function(x){ y = cor.test(as.numeric(eset_Tumor[probe, ]), as.numeric(x), alternative = "two.sided", method = met);return(y$p.value)})
            cor$probe_ID <- rownames(cor)
            cor_0.3 <- cor[abs(cor[, GSE_cor]) > 0.3 & cor[, GSE_p] < 0.05, ]
            cor_0.3 <- na.omit(cor_0.3)
            cor_SYMBOL <- as.data.frame(getSYMBOL(rownames(cor_0.3), "hgu133plus2.db"))
            cor_SYMBOL$probe_ID <- rownames(cor_SYMBOL)
            cor_0.3 <- merge(cor_0.3, cor_SYMBOL, by.x = "probe_ID")
            colnames(cor_0.3)[4] <- "SYMBOL"
            #cor_0.3 <- cor_0.3[order(cor_0.3[, 2], decreasing = T), ]
            return(cor_0.3)
}

get_cor_3.0 <- function(eset, group, T, met, probe1, probe2){eset_Tumor <- eset[, grep(T, group)] # "T" or "Tumor"
            cor <- cor(as.numeric(eset_Tumor[probe1, ]), as.numeric(eset_Tumor[probe2, ]), method = met, use = "pair")
            # "pearson" (default), "kendall", or "spearman"
            cor_test <- cor.test(as.numeric(eset_Tumor[probe1, ]), as.numeric(eset_Tumor[probe2, ]), alternative = "two.sided", method = met)
            cor_0.3 <- list(cor, cor_test$p.value)
            return(cor_0.3)
}




system.time(GSE45267_OLM_cor <- get_cor2.0(GSE45267_expression, GSE45267_group, "232230_at"))



#################################################### GSE45267 dif plot ############################################

#get_plot <- function(probe_name, gene_name){
#  GSE_data <- get_SYMBOL_data(probe_name, gene_name, GSE45267_expression, GSE45267_group)
#  GSE_plot <- ggplot(GSE_data, aes(y = UC002TVK, x = group , colour = group, shape = group))+
#  geom_boxplot(width = 0.5)+scale_x_discrete(labels = c("Tumor", "Normal"), limits = c("T", "N"))+
#  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
#  theme_linedraw()+labs(title = paste("GSE45267", gene_name, sep = "-"), x= "", y="Expression (log2)")+
#  annotate("text",x=0.5,y=5.5,label = paste("Fold_change =", signif(GSE45267_dif[probe_name, 1],4), sep = ""),hjust=0)+
#  annotate("text", x=0.5, y=5.2, label = paste("P.Val =", signif(GSE45267_dif[probe_name, 4],4), sep = ""),hjust=0)
#  return(GSE_plot)
#}

#GSE45267_ZEB1_plot <- get_plot("229090_at","ZEB1", GSE45267_expression, group, "GSE45267", GSE45267_dif)
#GSE45267_PCAT6_plot <- get_plot("224945_at", "ENST00000417262", GSE45267_expression, group, "GSE45267", GSE45267_dif)
#GSE45267_IGLL5_plot <- get_plot("215121_x_at", "NR_033661", GSE45267_expression, group, "GSE45267", GSE45267_dif)
#GSE45267_AK129699_plot <- get_plot("238619_at", "AK129699")
#GSE45267_BC043009_plot <- get_plot("230930_at", "BC043009")
#GSE45267_ENST00000449772 <- get_plot("1563298_at", "ENST00000449772")
#GSE45267_UC002TVK <- get_plot("210662_at", "UC002TVK")
################################################## GSE55092 dif plot ###########

#setwd("/home/yhy/liver_cancer/GSE55092/")
#cel.files <- list.files(pattern = ".cel$", ignore.case = T, full.names = T)
#GSE55092_rawdata <- read.affybatch(cel.files)
#GSE55092_gcrma <- gcrma::gcrma(GSE55092_rawdata)
#GSE55092_expression <- exprs(GSE55092_gcrma)
#GSE55092_expression <- as.data.frame(GSE55092_expression)



GSE55092_group <- apply(as.array(colnames(GSE55092_expression)), 1, function(x)substr(x, 12, nchar(x)-4))

GSE55092_dif <- get_flodchange(GSE55092_group, GSE55092_expression, 1)
#rownames(MCU_expression) %in% rownames(GSE55092_dif)
#rownames(OLMALINC_expression) %in% rownames(GSE55092_dif)
#GSE55092_dif[rownames(SCD_expression),]


library(ggplot2)
library(reshape2)


#get_plot <- function(probe_name, gene_name){
#  GSE_data <- get_SYMBOL_data(probe_name, GSE55092_expression, GSE55092_group)
#  GSE_plot <- ggplot(GSE_data, aes(y = SYMBOL, x = group , colour = group, shape = group))+
#  geom_boxplot(width = 0.5)+scale_x_discrete(limits = c("Tumor", "Normal"))+
#  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
#  theme_linedraw()+labs(title = paste("GSE55092", gene_name, sep = "-"), x= "", y="Expression (log2)")+
#  annotate("text",x=0.5,y=5.5,label = paste("Fold_change =", signif(GSE55092_dif[probe_name, 1],4), sep = ""),hjust=0)+
#  annotate("text", x=0.5, y=5.2, label = paste("P.Val =", signif(GSE55092_dif[probe_name, 4],4), sep = ""),hjust=0)
#  return(GSE_plot)
#}

#GSE55092_PCAT6_plot <- get_plot("224945_at", "ENST00000417262", GSE55092_expression, GSE55092_group, "GSE55092", GSE55092_dif)
#GSE55092_ZEB1_plot <- get_plot("229090_at","ZEB1", GSE55092_expression, GSE55092_group, "GSE55092", GSE55092_dif)
#GSE55092_IGLL5_plot <- get_plot("215121_x_at", "NR_033661", GSE55092_expression, GSE55092_group, "GSE55092", GSE55092_dif)
#GSE55092_AK129699_plot <- get_plot("238619_at", "AK129699")
#GSE55092_BC043009_plot <- get_plot("230930_at", "BC043009")
#GSE55092_ENST00000449772 <- get_plot("1563298_at", "ENST00000449772")


#for (i in c("GSE45267_PCAT6_plot", "GSE45267_ZEB1_plot", "GSE45267_IGLL5_plot", "GSE45267_AK129699_plot", "GSE55092_PCAT6_plot", "GSE55092_ZEB1_plot", "GSE55092_IGLL5_plot", "GSE55092_AK129699_plot", "GSE45267_BC043009_plot","GSE45267_ENST00000449772", "GSE55092_BC043009_plot", "GSE55092_ENST00000449772")) {
#ggsave(get(i), filename = paste("~/liver_cancer/test_plot/", i, ".jpeg", sep = ""), device = "jpeg")  
  
#}


#################################### cor ###########################

GSE45267_cor_pearson <- get_cor("GSE45267", GSE45267_expression, group, "T", "pearson") #2175
GSE45267_cor_spearman <- get_cor("GSE45267", GSE45267_expression, group, "T", "spearman")#2788
GSE45267_cor_kendall <- get_cor("GSE45267", GSE45267_expression, group, "T", "kendall")#360
#GSE45267_cor_SYMBOL <- Reduce(intersect, list(GSE45267_cor_pearson$SYMBOL, GSE45267_cor_spearman$SYMBOL, GSE45267_cor_kendall$SYMBOL))


GSE55092_cor_pearson <- get_cor("GSE55092", GSE55092_expression, GSE55092_group, "Tumor", "pearson")#4610
GSE55092_cor_spearman <- get_cor("GSE55092", GSE55092_expression, GSE55092_group, "Tumor", "spearman")#5524
GSE55092_cor_kendall <- get_cor("GSE55092", GSE55092_expression, GSE55092_group, "Tumor", "kendall")#1221
#GSE55092_cor_SYMBOL <- Reduce(intersect, list(GSE55092_cor_pearson$SYMBOL, GSE55092_cor_spearman$SYMBOL, GSE55092_cor_kendall$SYMBOL))


GSE45436_cor_pearson <- get_cor("GSE45436", GSE45436_expression, GSE45436_group, "T", "pearson")#1626
GSE45436_cor_spearman <- get_cor("GSE45436", GSE45436_expression, GSE45436_group, "T", "spearman")#1427
GSE45436_cor_kendall <- get_cor("GSE45436", GSE45436_expression, GSE45436_group, "T", "kendall")#75
#GSE45436_cor_SYMBOL <- Reduce(intersect, list(GSE45436_cor_pearson$SYMBOL, GSE45436_cor_spearman$SYMBOL, GSE45436_cor_kendall$SYMBOL))


cor_pearson <- Reduce(intersect, list(GSE45436_cor_pearson$SYMBOL, GSE55092_cor_pearson$SYMBOL, GSE45267_cor_pearson$SYMBOL)) #203
cor_spearman <- Reduce(intersect, list(GSE45267_cor_spearman$SYMBOL, GSE55092_cor_spearman$SYMBOL, GSE45436_cor_spearman$SYMBOL)) #213
cor_kendall <- Reduce(intersect, list(GSE45267_cor_kendall$SYMBOL, GSE55092_cor_kendall$SYMBOL, GSE45436_cor_kendall$SYMBOL))

something_cor <- intersect(cor_pearson, cor_spearman)

cor_pearson_probe <- Reduce(intersect, list(GSE45436_cor_pearson$probe_ID, GSE45267_cor_pearson$probe_ID)) #919
cor_spearman_probe <- Reduce(intersect, list(GSE45267_cor_spearman$probe_ID, GSE45436_cor_spearman$probe_ID)) #779
#cor_kendall_probe <- Reduce(intersect, list(GSE45267_cor_kendall$probe_ID, GSE55092_cor_kendall$probe_ID, GSE45436_cor_kendall$probe_ID))




pearson_probe <- Reduce(function(x, y) merge(x, y, by = "probe_ID",), list(GSE45436_cor_pearson[GSE45436_cor_pearson$probe_ID %in% cor_pearson_probe, ], GSE45267_cor_pearson[GSE45267_cor_pearson$probe_ID %in% cor_pearson_probe, ]))
spearman_probe <- Reduce(function(x, y) merge(x, y, by = "probe_ID",), list(GSE45436_cor_spearman[GSE45436_cor_spearman$probe_ID %in% cor_spearman_probe, ], GSE45267_cor_spearman[GSE45267_cor_spearman$probe_ID %in% cor_spearman_probe, ]))


pearson_comm_up <- pearson_probe[pearson_probe$SYMBOL.x %in% intersect(dif$SYMBOL[dif$logFC > 0], pearson_probe$SYMBOL.x), ]
pearson_comm_down <- pearson_probe[pearson_probe$SYMBOL.x %in% intersect(dif$SYMBOL[dif$logFC < 0], pearson_probe$SYMBOL.x), ]
spearman_comm_up <- spearman_probe[spearman_probe$SYMBOL.x %in% intersect(dif$SYMBOL[dif$logFC > 0], spearman_probe$SYMBOL.x), ]
spearman_comm_down <- spearman_probe[spearman_probe$SYMBOL.x %in% intersect(dif$SYMBOL[dif$logFC < 0], spearman_probe$SYMBOL.x), ]

#HNRNPK  "200097_s_at", "200775_s_at"
#MTHFD1  "1561125_at"  "202309_at"   "225518_x_at" "239846_at"

#ZEB1 "210875_s_at", "212758_s_at", "212764_at", "239952_at" 
#ZEB2 "203603_s_at", "228333_at", "233031_at", "233033_at", "235593_at"


GSE45267_expression_HNRNPK <- GSE45267_expression[c("200097_s_at", "200775_s_at"), ]
GSE45436_expression_HNRNPK <- GSE45436_expression[c("200097_s_at", "200775_s_at"), ]

GSE45267_cor_pearson_HNRNPK_1 <- get_cor_2.0("GSE45267", GSE45267_expression, group, "T", "pearson", "200097_s_at") 
GSE45267_cor_spearman_HNRNPK_1 <- get_cor_2.0("GSE45267", GSE45267_expression, group, "T", "spearman", "200097_s_at")
GSE45267_cor_pearson_HNRNPK_2 <- get_cor_2.0("GSE45267", GSE45267_expression, group, "T", "pearson", "200775_s_at") 
GSE45267_cor_spearman_HNRNPK_2 <- get_cor_2.0("GSE45267", GSE45267_expression, group, "T", "spearman", "200775_s_at")

GSE45436_cor_pearson_HNRNPK_1 <- get_cor_2.0("GSE45436", GSE45436_expression, group, "T", "pearson", "200097_s_at") 
GSE45436_cor_spearman_HNRNPK_1 <- get_cor_2.0("GSE45436", GSE45436_expression, group, "T", "spearman", "200097_s_at")
GSE45436_cor_pearson_HNRNPK_2 <- get_cor_2.0("GSE45436", GSE45436_expression, group, "T", "pearson", "200775_s_at") 
GSE45436_cor_spearman_HNRNPK_2 <- get_cor_2.0("GSE45436", GSE45436_expression, group, "T", "spearman", "200775_s_at")

HNRPK_pearson_probe <- Reduce(intersect, list(GSE45267_cor_pearson_HNRNPK_1$probe_ID, GSE45267_cor_pearson_HNRNPK_2$probe_ID, GSE45436_cor_pearson_HNRNPK_1$probe_ID, GSE45436_cor_pearson_HNRNPK_2$probe_ID))
HNRPK_spearman_probe <- Reduce(intersect, list(GSE45267_cor_spearman_HNRNPK_1$probe_ID, GSE45267_cor_spearman_HNRNPK_2$probe_ID, GSE45436_cor_spearman_HNRNPK_1$probe_ID, GSE45436_cor_spearman_HNRNPK_2$probe_ID))

HNRNPK_pearson <- Reduce(function(x, y) merge(x, y, by = "probe_ID",), list(GSE45267_cor_pearson_HNRNPK_1[GSE45267_cor_pearson_HNRNPK_1$probe_ID %in% HNRPK_pearson_probe, ], GSE45267_cor_pearson_HNRNPK_2[GSE45267_cor_pearson_HNRNPK_2$probe_ID %in% HNRPK_pearson_probe, ]
                                                                            , GSE45436_cor_pearson_HNRNPK_1[GSE45436_cor_pearson_HNRNPK_1$probe_ID %in% HNRPK_pearson_probe, ], GSE45436_cor_pearson_HNRNPK_2[GSE45436_cor_pearson_HNRNPK_2$probe_ID %in% HNRPK_pearson_probe, ]))

HNRNPK_spearman <- Reduce(function(x, y) merge(x, y, by = "probe_ID",), list(GSE45267_cor_spearman_HNRNPK_1[GSE45267_cor_spearman_HNRNPK_1$probe_ID %in% HNRPK_spearman_probe, ], GSE45267_cor_spearman_HNRNPK_2[GSE45267_cor_spearman_HNRNPK_2$probe_ID %in% HNRPK_spearman_probe, ]
                                                                            , GSE45436_cor_spearman_HNRNPK_1[GSE45436_cor_spearman_HNRNPK_1$probe_ID %in% HNRPK_spearman_probe, ], GSE45436_cor_spearman_HNRNPK_2[GSE45436_cor_spearman_HNRNPK_2$probe_ID %in% HNRPK_spearman_probe, ]))








#########################################
library(GOstats)

hgu133plus2_EG <- as.character(unique(na.omit(getEG(rownames(GSE45436_expression), "hgu133plus2.db"))))

GO_analyse <- function(co_expression_set){
            pearson_EG <- as.character(unique(na.omit(getEG(co_expression_set$probe_ID, "hgu133plus2.db"))))
            params <- new("GOHyperGParams", geneIds = pearson_EG, universeGeneIds = hgu133plus2_EG, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.001, conditional = F, testDirection = "over")
            hgOver <- hyperGTest(params)
            bp <- summary(hgOver)
            return(bp)
            }

pearson_GO <- GO_analyse(pearson_probe)
pearson_UP_GO <- GO_analyse(pearson_probe[pearson_probe$GSE45436_cor > 0 & pearson_probe$GSE45267_cor > 0, ])
pearson_DOWN_GO <- GO_analyse(pearson_probe[pearson_probe$GSE45436_cor < 0 & pearson_probe$GSE45267_cor < 0, ])


spearman_GO <- GO_analyse(spearman_probe)
spearman_UP_GO <- GO_analyse(spearman_probe[spearman_probe$GSE45436_cor > 0 & spearman_probe$GSE45267_cor > 0, ])
spearman_DOWN_GO <- GO_analyse(spearman_probe[spearman_probe$GSE45436_cor < 0 & spearman_probe$GSE45267_cor < 0, ])

for (i in c("pearson_GO", "pearson_UP_GO", "pearson_DOWN_GO", "spearman_GO", "spearman_UP_GO", "spearman_DOWN_GO")) {
  write.csv(get(i), paste("/home/yhy/liver_cancer/", "OLM_", i, ".csv", sep = ""))
}


HNRNPK_pearson_GO <- GO_analyse(HNRNPK_pearson)
HNRNPK_pearson_UP_GO <- GO_analyse(HNRNPK_pearson[HNRNPK_pearson$GSE45267_cor.x > 0 & HNRNPK_pearson$GSE45267_cor.y > 0 & HNRNPK_pearson$GSE45436_cor.x > 0 & HNRNPK_pearson$GSE45436_cor.y > 0, ])
HNRNPK_pearson_DOWN_GO <- GO_analyse(HNRNPK_pearson[HNRNPK_pearson$GSE45267_cor.x < 0 & HNRNPK_pearson$GSE45267_cor.y < 0 & HNRNPK_pearson$GSE45436_cor.x < 0 & HNRNPK_pearson$GSE45436_cor.y < 0, ])

HNRNPK_spearman_GO <- GO_analyse(HNRNPK_spearman)
HNRNPK_spearman_UP_GO <- GO_analyse(HNRNPK_spearman[HNRNPK_spearman$GSE45267_cor.x > 0 & HNRNPK_spearman$GSE45267_cor.y > 0 & HNRNPK_spearman$GSE45436_cor.x > 0 & HNRNPK_spearman$GSE45436_cor.y > 0, ])
HNRNPK_spearman_DOWN_GO <- GO_analyse(HNRNPK_spearman[HNRNPK_spearman$GSE45267_cor.x < 0 & HNRNPK_spearman$GSE45267_cor.y < 0 & HNRNPK_spearman$GSE45436_cor.x < 0 & HNRNPK_spearman$GSE45436_cor.y < 0, ])

for (i in c("HNRNPK_pearson_GO", "HNRNPK_pearson_UP_GO", "HNRNPK_pearson_DOWN_GO", "HNRNPK_spearman_GO", "HNRNPK_spearman_UP_GO", "HNRNPK_spearman_DOWN_GO")) {
  write.csv(get(i), paste("/home/yhy/liver_cancer/", i, ".csv", sep = ""))
}




GO_analyse_cor <- function(pearson_probe, spearman_probe, GSE_cor){
                      pearson_GO <- GO_analyse(pearson_probe)
                      pearson_UP_GO <- GO_analyse(pearson_probe[pearson_probe[, GSE_cor] > 0 , ])
                      pearson_DOWN_GO <- GO_analyse(pearson_probe[pearson_probe[, GSE_cor] < 0 , ])
                      
                      spearman_GO <- GO_analyse(spearman_probe)
                      spearman_UP_GO <- GO_analyse(spearman_probe[spearman_probe[, GSE_cor] > 0 , ])
                      spearman_DOWN_GO <- GO_analyse(spearman_probe[spearman_probe[, GSE_cor] < 0 , ])
                      return(list(pearson_GO, pearson_UP_GO, pearson_DOWN_GO, spearman_GO, spearman_UP_GO, spearman_DOWN_GO))
}

GSE45267_cor_GO <- GO_analyse_cor(GSE45267_cor_pearson, GSE45267_cor_spearman, "GSE45267_cor")
GSE45436_cor_GO <- GO_analyse_cor(GSE45436_cor_pearson, GSE45436_cor_spearman, "GSE45436_cor")

GSE_cor_GO_name <- c("pearson_GO", "pearson_UP_GO", "pearson_DOWN_GO", "spearman_GO", "spearman_UP_GO", "spearman_DOWN_GO")
for (i in 1:length(GSE_cor_GO_name)) {
  write.csv(GSE45436_cor_GO[[i]], paste("/home/yhy/liver_cancer/", "GSE45436_", GSE_cor_GO_name[i], ".csv", sep = ""))
  
}

for (i in 1:length(GSE_cor_GO_name)) {
  write.csv(GSE45267_cor_GO[[i]], paste("/home/yhy/liver_cancer/", "GSE45267_", GSE_cor_GO_name[i], ".csv", sep = ""))
  
}

######################################################

GO_analyse_dif <- function(GSE_dif){
  EG_Total <- as.character(unique(na.omit(getEG(rownames(GSE_dif), "hgu133plus2.db"))))
  EG_UP <- as.character(unique(na.omit(getEG(rownames(GSE_dif[GSE_dif$logFC > 0, ]), "hgu133plus2.db"))))
  EG_DOWN <- as.character(unique(na.omit(getEG(rownames(GSE_dif[GSE_dif$logFC < 0, ]), "hgu133plus2.db"))))
  params_Total <- new("GOHyperGParams", geneIds = EG_Total, universeGeneIds = hgu133plus2_EG, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.001, conditional = F, testDirection = "over")
  params_UP <- new("GOHyperGParams", geneIds = EG_UP, universeGeneIds = hgu133plus2_EG, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.001, conditional = F, testDirection = "over")
  params_DOWN <- new("GOHyperGParams", geneIds = EG_DOWN, universeGeneIds = hgu133plus2_EG, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.001, conditional = F, testDirection = "over")
  hgOver_Total <- hyperGTest(params_Total)
  hgOver_UP <- hyperGTest(params_UP)
  hgOver_DOWN <- hyperGTest(params_DOWN)
  bp_Total <- summary(hgOver_Total)
  bp_UP <- summary(hgOver_UP)
  bp_DOWN <- summary(hgOver_DOWN)
  return(list(bp_Total, bp_UP, bp_DOWN))
}

GSE45436_dif_GO <- GO_analyse_dif(GSE45436_dif)
GSE45436_dif_GO_Total <- GSE45436_dif_GO[[1]]
GSE45436_dif_GO_UP <- GSE45436_dif_GO[[2]]
GSE45436_dif_GO_DOWN <- GSE45436_dif_GO[[3]]

GSE45267_dif_GO <- GO_analyse_dif(GSE45267_dif)
GSE45267_dif_GO_Total <- GSE45267_dif_GO[[1]]
GSE45267_dif_GO_UP <- GSE45267_dif_GO[[2]]
GSE45267_dif_GO_DOWN <- GSE45267_dif_GO[[3]]

for (i in c("GSE45436_dif_GO_Total", "GSE45436_dif_GO_UP", "GSE45436_dif_GO_DOWN", "GSE45267_dif_GO_Total", "GSE45267_dif_GO_UP", "GSE45267_dif_GO_DOWN"))
{
  file_save <- get(i)
  write.csv(file_save, paste("/home/yhy/liver_cancer/", i, ".csv", sep = ""))
}  


##########################################################
#ZEB1 "210875_s_at", "212758_s_at", "212764_at", "239952_at" 
#ZEB2 "203603_s_at", "228333_at", "233031_at", "233033_at", "235593_at"

ZBE_probe <- c("210875_s_at", "212758_s_at", "212764_at", "239952_at", "203603_s_at", "228333_at", "233031_at", "233033_at", "235593_at")
GSE45436_ZEB <- GSE45436_expression[ZBE_probe, ]
apply(GSE45436_ZEB, 1, max)
apply(GSE45436_ZEB, 1, mean)
#ZEB1 "212764_at" , ZEB2 "203603_s_at"

GSE45436_ZEB1_cor <- get_cor_3.0(GSE45436_expression, GSE45436_group, "T", "spearman", "1568611_at", "212764_at")
GSE45436_ZEB2_cor <- get_cor_3.0(GSE45436_expression, GSE45436_group, "T", "spearman", "1568611_at", "203603_s_at")
GSE45267_ZEB1_cor <- get_cor_3.0(GSE45267_expression, group, "T", "spearman", "1568611_at", "212764_at")
GSE45267_ZEB2_cor <- get_cor_3.0(GSE45267_expression, group, "T", "spearman", "1568611_at", "203603_s_at")

#P4HA2 "1555027_at", "202733_at"
GSE45436_ZEB1_cor <- get_cor_3.0(GSE45436_expression, GSE45436_group, "T", "spearman", "1568611_at", "1555027_at")
GSE45436_ZEB2_cor <- get_cor_3.0(GSE45436_expression, GSE45436_group, "T", "spearman", "1568611_at", "202733_at")###
GSE45267_ZEB1_cor <- get_cor_3.0(GSE45267_expression, group, "T", "spearman", "1568611_at", "1555027_at")
GSE45267_ZEB2_cor <- get_cor_3.0(GSE45267_expression, group, "T", "spearman", "1568611_at", "202733_at")####


GSE45436_P4HA2_data <- get_SYMBOL_data("202733_at", "P4HA2", GSE45436_expression, GSE45436_group)
ggplot(GSE45436_P4HA2_data, aes(x=group, y= P4HA2, colour = group, shape = group))+
  geom_boxplot(width = 0.5)+
  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
  theme_linedraw()+labs(title = "GSE45436-P4HA2", x= "", y="expression (log2)")+
  annotate("text",x=0.5,y=10.5,label = paste("Fold_change =", signif(GSE45436_dif["202733_at", 1],4), sep = ""),hjust=0)+
  annotate("text", x=0.5, y=10.2, label = paste("P.Val =", signif(GSE45436_dif["202733_at", 4],4), sep = ""),hjust=0)

GSE45267_P4HA2_data <- get_SYMBOL_data("202733_at", "P4HA2", GSE45267_expression, group)
GSE45267_P4HA2_plot <- ggplot(GSE45267_P4HA2_data, aes(x=group, y=P4HA2 , colour = group, shape = group))+
  geom_boxplot(width = 0.5)+
  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
  theme_linedraw()+labs(title = "GSE45267-P4HA2", x= "", y="expression (log2)")+
  annotate("text",x=0.5,y=9.5,label = paste("Fold_change =", signif(GSE45267_dif["202733_at", 1],4), sep = ""),hjust=0)+
  annotate("text", x=0.5, y=9.2, label = paste("P.Val =", signif(GSE45267_dif["202733_at", 4],4), sep = ""),hjust=0)

save(GSE45267_BC_data, GSE55092_BC_data, GSE45436_BC_data, GSE45267_dif, GSE55092_dif, GSE45436_dif, file = "/home/yhy/liver_cancer/BC_plotdata.RData")
#########################################################

#survival

#########################################################
library(survival)
library(oligo)
library(affy)
GPL10558 <- read.delim("~/liver_cancer/GPL10558/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt", header = T, stringsAsFactors = F)
GSE40873_cel <- list.celfiles("~/liver_cancer/GSE40873/")
setwd("~/liver_cancer/GSE40873/")
GSE40873_raw <- read.affybatch(GSE40873_cel)
GSE40873_gcrma <- gcrma::gcrma(GSE40873_raw)
GSE40873_exprs <- exprs(GSE40873_gcrma)
GSE40873_clinical <- read.delim("~/liver_cancer/GSE40873/GSE40873_series_matrix.txt", header = T, stringsAsFactors = F, sep = ",")#11 rows is survial 
rownames(GSE40873_clinical)[c(11,40:nrow(GSE40873_clinical))] <- c("survival", GSE40873_clinical$X.Sample_title[40:nrow(GSE40873_clinical)])
GSE40873_clinical <- GSE40873_clinical[,-1]
GSE40873_clinical["survival_day",] <- lapply(GSE40873_clinical["survival",], FUN = function(x){strsplit(x, ":")[[1]][2]})
GSE40873_clinical["events",] <- lapply(GSE40873_clinical[12,], FUN = function(x){strsplit(x, ":")[[1]][2]})
GSE40873_clinical["MO_status",] <- grepl("^multicentric occurrence", GSE40873_clinical["events",])

LincRNA_survival <- t(GSE40873_clinical[c("survival_day",probe_gene$probe),])         #"232230_at", "229090_at", "224945_at", "215121_x_at", "238619_at", "230930_at", "1563298_at", "239503_at"),])
LincRNA_survival <- as.data.frame(apply(LincRNA_survival, 2, FUN = function(x){as.numeric(x)}))
LincRNA_survival$MO_status <- t(GSE40873_clinical["MO_status",])
LincRNA_survival$MO_status <- as.logical(t(GSE40873_clinical["MO_status",]))


library(survminer)
probe_gene <- data.frame(list(probe = c("229090_at", "224945_at", "215121_x_at", "238619_at", "230930_at", "1563298_at", "239503_at"), gene = c("ZEB1_AS1", "ENST00000417262", "NR_033661", "AK129699", "BC043009", "ENST00000449772", "BC041955")), stringsAsFactors = F)
probe_gene[8, ] <- c("232230_at", "OLMALINC")
probe_gene[9:10,] <- c("200831_s_at", "200832_s_at", "SCD1", "SCD2")
probe_gene[11:17,] <- c("210795_s_at", "217797_at", "217494_s_at", "239153_at", "212495_at", "210662_at", "1561542_at", "MEG3", "UFC1", "PTENP1", "HOTAIR", "BC015378", "UC002TVK", "AK056988")
probe_gene[18:19,] <- c("1568611_at", "225457_s_at", "BC013423", "OLMALINC")




for (i in 1:nrow(probe_gene)) {
      y <- Surv(LincRNA_survival$survival_day, LincRNA_survival$MO_status)
      probe_cutoff <- cutoff.survival(LincRNA_survival, probe_gene$probe[i])
      group <- LincRNA_survival[, probe_gene$probe[i]] <= probe_cutoff$cutoff
      kmfit <- survfit(y~group)
      jpeg(paste("~/liver_cancer/GSE40873_survival/GSE40873-", probe_gene$gene[i], ".jpeg", sep = ""))
      print(ggsurvplot(kmfit, pval = T, legend.title = paste("GSE40873-", probe_gene$gene[i], sep = ""), legend.labs = c("high_expression", "low_expression"), risk.table = T))
      dev.off()
      print(i)
}

cor(as.numeric(LincRNA_survival$`232230_at`), as.numeric(LincRNA_survival$`200831_s_at`), method = "pearson")
cor(as.numeric(LincRNA_survival$`225457_s_at`), as.numeric(LincRNA_survival$`200831_s_at`), method = "pearson")
cor(as.numeric(LincRNA_survival$`232230_at`), as.numeric(LincRNA_survival$`200832_s_at`), method = "pearson")
cor(as.numeric(LincRNA_survival$`225457_s_at`), as.numeric(LincRNA_survival$`200832_s_at`), method = "pearson")


#myplot <- ggsurvplot(kmfit, pval = T, legend.title = paste("GSE40873-", probe_gene$gene[i], sep = ""), legend.labs = c("high_expression", "low_expression"), risk.table = T)


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#Test for survival

###########################################























##########################################################

#   QC   2016-12-7############

##########################################################

library(affyPLM)
library(simpleaffy)
GSE45436_rawdata.qc <- qc(GSE45436_rawdata)
GSE45436_Pset <- fitPLM(GSE45436_rawdata)

for (i in 1:ncol(GSE45436_expression)){
  jpeg(filename = paste("~/liver_cancer/GSE45436/resids/", colnames(GSE45436_expression)[i], sep = ""))
  image(GSE45436_Pset, type = "resids", which = i)
  dev.off()
}

jpeg(filename = "~/liver_cancer/GSE45436/resids/QC.jpeg", width = 4800, height = 20000, pointsize = 120)
plot(GSE45436_rawdata.qc)
dev.off()
jpeg(filename = "~/liver_cancer/GSE45436/resids/RLE.jpeg", width = 4800, height = 2000, pointsize = 120)
Mbox(GSE45436_Pset, main = "RLE", ylim = c(-1,1))
dev.off()
jpeg(filename = "~/liver_cancer/GSE45436/resids/NUSE.jpeg", width = 4800, height = 2000, pointsize = 120)
boxplot(GSE45436_Pset, main = "NUSE", ylim = c(0.96, 1.06))
dev.off()

GSE45436_qc.probes <- as.data.frame(GSE45436_rawdata.qc@qc.probes)
#GSE45436_degradate <- rownames(GSE45436_qc.probes)[GSE45436_qc.probes$`AFFX-HSAC07/X00351_3_at`-GSE45436_qc.probes$`AFFX-HSAC07/X00351_5_at` > 3 | GSE45436_qc.probes$`AFFX-HUMGAPDH/M33197_3_at`-GSE45436_qc.probes$`AFFX-HUMGAPDH/M33197_5_at` > 1.25]
GSE45436_degradate <- rownames(GSE45436_qc.probes)[GSE45436_qc.probes$`AFFX-HSAC07/X00351_3_at`-GSE45436_qc.probes$`AFFX-HSAC07/X00351_5_at` > 3] 


#GSE45436_NUSE <- as.data.frame(NUSE(GSE45436_Pset, type = "stats"))
#GSE45436_NUSEhigh <- colnames(GSE45436_NUSE)[GSE45436_NUSE[1,] > 1.02]
#intersect(GSE45436_degradate, GSE45436_NUSEhigh)
#GSE45436_delet <- union(GSE45436_degradate, GSE45436_NUSEhigh)



GSE45436_deg_affy <- AffyRNAdeg(GSE45436_rawdata)
plotAffyRNAdeg(GSE45436_deg_affy)
#save.image("~/liver_cancer/microarray.RData")
#load("~/liver_cancer/microarray.RData")


GSE45436_delet <- union(list.celfiles("~/liver_cancer/GSE45436/GSE45434_repeat/"), GSE45436_degradate)
GSE45436_celfile <- list.celfiles("~/liver_cancer/GSE45436/")
setwd("~/liver_cancer/GSE45436/")
GSE45436_celfile <- setdiff(GSE45436_celfile, GSE45436_delet)
GSE45436_rawdata <- read.affybatch(GSE45436_celfile) 
GSE45436_new.qc <- qc(GSE45436_rawdata)
GSE45436_new.Pset <- fitPLM(GSE45436_rawdata)
GSE45436_new_NUSE <- as.data.frame(NUSE(GSE45436_new.Pset, type = "stats"))
GSE45436_new_NUSE.high <- colnames(GSE45436_new_NUSE)[GSE45436_new_NUSE[1,] > 1.02]



GSE45436_gcrma <- gcrma::gcrma(GSE45436_rawdata)
GSE45436_new_expression <- exprs(GSE45436_gcrma)

GSE45436_color <- data.frame(color = rep("blue", length(colnames(GSE45436_new_expression))), stringsAsFactors = F)
rownames(GSE45436_color) <- colnames(GSE45436_new_expression)
GSE45267_celfile <- list.celfiles("~/liver_cancer/GSE45267/")
GSE45435_celfile <- setdiff(GSE45436_celfile, GSE45267_celfile)
GSE45436_color[GSE45435_celfile,1] <- "green"
GSE45436_color[GSE45436_new_NUSE.high,] <- "red"
GSE45436_qc_hc <- get_qc_hc(GSE45436_new_expression)
library(dendextend)
GSE45436_hc_forplot <- as.dendrogram(GSE45436_qc_hc)
labels_colors(GSE45436_hc_forplot) <- GSE45436_color$color[order.dendrogram(GSE45436_hc_forplot)]
par(mai = c(2,0.1,0.1,0.1), fin = c(10, 6.99999))
jpeg("~/liver_cancer/GSE45436/GSE45436_qchc.jpeg", width = 1000, height = 1500)
plot(GSE45436_hc_forplot)
dev.off()


library(estimate)
library(hgu133plus2.db)
GSE45436_symbol_expression <- data.frame(probe_symbol, GSE45436_new_expression[rownames(probe_symbol),])
GSE45436_symbol_expression <- aggregate(.~SYMBOL, data = GSE45436_symbol_expression, max)
rownames(GSE45436_symbol_expression) <- GSE45436_symbol_expression$SYMBOL
GSE45436_symbol_expression <- GSE45436_symbol_expression[,-1]
#OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package="estimate")
write.table(GSE45436_symbol_expression, "~/liver_cancer/GSE45436/GSE45436_estimate_input.txt", sep = "\t", quote = F)
filterCommonGenes(input.f = "~/liver_cancer/GSE45436/GSE45436_estimate_input.txt", output.f = "~/liver_cancer/GSE45436/GSE45436_estimate_output.gct", id = "GeneSymbol")
estimateScore("~/liver_cancer/GSE45436/GSE45436_estimate_output.gct", "~/liver_cancer/GSE45436/GSE45436_estimate_score.gct", platform="affymetrix")
plotPurity(scores = "~/liver_cancer/GSE45436/GSE45436_estimate_score.gct", platform = "affymetrix")

###################################################
####delete NUSE.high and estimate tumor low########
###################################################
GSE45436_celfile <- setdiff(GSE45436_celfile, GSE45436_new_NUSE.high)
estimate_delet_name <- c("15T", "377T", "89T", "HCC_30T", "1717N", "001T.2TT")
GSE45436_celfile <- setdiff(GSE45436_celfile, GSE45436_celfile[unlist(lapply(estimate_delet_name, FUN = function(x) grep(x, GSE45436_celfile)))])
GSE45436_rawdata <- read.affybatch(GSE45436_celfile)
GSE45436_gcrma <- gcrma::gcrma(GSE45436_rawdata)
GSE45436_final_expression <- exprs(GSE45436_gcrma)
GSE45436_final_expression <- as.data.frame(GSE45436_final_expression)
GSE45436_final_hc <- get_qc_hc(GSE45436_final_expression)
plot(GSE45436_final_hc)

par(mai = c(0.1,0.1,0.1,0.1), fin = c(7, 5.99999))


GSE45436_final_group <- unlist(lapply(colnames(GSE45436_final_expression), FUN = function(x)str_right(x,5)))
GSE45436_final_group <- sub("T", "Tumor", GSE45436_final_group)
GSE45436_final_group <- sub("N", "Normal", GSE45436_final_group)
GSE45436_final_dif <- get_flodchange(GSE45436_final_group, GSE45436_final_expression, 1)





for(i in 1:nrow(probe_gene)){
  plot_name <- paste(probe_gene$gene[i], "plot", sep = "_")
    assign(plot_name, get_plot(probe_gene$probe[i], probe_gene$gene[i], GSE45436_final_expression, GSE45436_final_group, GSE45436_final_dif, "GSE45436")) 
  ggsave(paste("~/liver_cancer/GSE45436/GSE45436_final_dif/", plot_name, ".jpeg", sep = ""), plot = get(plot_name), device = "jpeg")
}

GSE45436_group <- sub("T", "Tumor", GSE45436_group)
GSE45436_group <- sub("N", "Normal", GSE45436_group)
for(i in 1:nrow(probe_gene)){
  plot_name <- paste(probe_gene$gene[i], "plot", sep = "_")
  assign(plot_name, get_plot(probe_gene$probe[i], probe_gene$gene[i], GSE45436_expression, GSE45436_group, GSE45436_dif, "GSE45436_old")) 
  ggsave(paste("~/liver_cancer/GSE45436/GSE45436_old_dif/", plot_name, ".jpeg", sep = ""), plot = get(plot_name), device = "jpeg")
}

make_script2.0(GSE45436_expression, GSE45436_group, "GSE45436_old")
make_script2.0(GSE45436_final_expression, GSE45436_final_group, "GSE45436_final")



get_cor_5.0 <- function(eset_name, eset, group, type, probe, met = "pearson"){
            eset <- eset[, grep(type, group)] # "T" or "Tumor"
            cor <- as.data.frame(apply(eset, 1, FUN = function(x){cor(as.numeric(eset[probe, ]), as.numeric(x), method = met,use = "pair")}))
            GSE_cor <- paste(eset_name, "cor", sep = "_")
            colnames(cor) <- GSE_cor
            GSE_p <- paste(eset_name, "P.value", sep = "_")
            cor[, GSE_p] <- apply(eset, 1, FUN = function(x){ y = cor.test(as.numeric(eset[probe, ]), as.numeric(x), alternative = "two.sided", method = met);return(y$p.value)})
            cor <- cor[abs(cor[, GSE_cor]) > 0.3 & cor[, GSE_p] < 0.05, ]
            cor <- na.omit(cor)
            
            coexpression_probe <- rownames(cor)
            coexpression_probe <- intersect(coexpression_probe, rownames(probe_symbol))
            eset <- data.frame(probe_symbol[coexpression_probe, ], eset[coexpression_probe, ])
            colnames(eset)[1] <- "SYMBOL"
            eset <- aggregate(.~SYMBOL, data = eset, max)
            rownames(eset) <- eset$SYMBOL
            eset <- eset[, -1]
            return(eset)
}


GSE45436_AK129699_expression <- get_cor_5.0("GSE45436", GSE45436_final_expression, GSE45436_final_group, "Tumor", "238619_at")
expres_gct(GSE45436_AK129699_expression, "test")
test <- data.frame("238619_at", "AK129699")
colnames(test) <- c("probe", "gene")
phenotype_cls(GSE45436_AK129699_expression, "test", test)


test2 <- data.frame("232230_at", "OLMLINC", stringsAsFactors = F)
colnames(test2) <- c("probe", "gene")
test2[2:3,] <- c("200831_s_at", "200832_s_at", "SCD1", "SCD2")
make_script2.0(GSE45436_final_expression, group = GSE45436_final_group, GSE_name = "GSE45436_final", probe_gene = test2)
make_script2.0(GSE45436_expression, GSE45436_group, "GSE45436_old", test2)
##############################################################################
##########
##########GSE6222#####################
##########
##############################################################################

GSE6222_cel <- list.celfiles("~/liver_cancer/GSE6222/")
setwd("~/liver_cancer/GSE6222/")
GSE6222_raw <- read.affybatch(GSE6222_cel)
GSE6222_qc <- qc(GSE6222_raw)
GSE6222_degradate <- c("GSM143553.CEL", "GSM143552.CEL", "GSM143551.CEL", "GSM143550.CEL")
GSE6222_cel <- setdiff(GSE6222_cel, GSE6222_degradate)
GSE6222_raw <- read.affybatch(GSE6222_cel)
GSE6222_PLM <- fitPLM(GSE6222_raw)
GSE6222_NUSE <- as.data.frame(NUSE(GSE6222_PLM, type = "stats"))
GSE6222_high <- colnames(GSE6222_NUSE)[GSE6222_NUSE[1,] > 1.02 | GSE6222_NUSE[1,] < 0.98] 
GSE6222_cel <- setdiff(GSE6222_cel, GSE6222_high)
GSE6222_raw <- read.affybatch(GSE6222_cel)
GSE6222_gcrma <- gcrma(GSE6222_raw)
GSE6222_expression <- exprs(GSE6222_gcrma)
GSE6222_hc <- get_qc_hc(GSE6222_expression)

  #test <- cbind(GSE45436_final_expression, GSE6222_expression)
#test <- get_qc_hc(test) ####GSM143548, GSM143549 maybe repeat


#GSE6222_ENST, GSE6222_refseq, GSE6222_Entrezg, GSE6222_lincRNA Version 20.0.0

#Version 21.0.0

GSE6222_raw@cdfName <- "hgu133plus2hsrefseqcdf"
GSE6222_rma <- rma(GSE6222_raw)
GSE6222_refseq_21 <- exprs(GSE6222_rma)

GSE6222_raw@cdfName <- "hgu133plus2hsentrezgcdf"
GSE6222_rma <- rma(GSE6222_raw)
GSE6222_Entrezg_21 <- exprs(GSE6222_rma)



ind <- rowSums(test == "A")



GSE55092_spike <- spikein_probe(GSE55092_rawdata)
GSE45436_spike <- spikein_probe(GSE45436_rawdata)
GSE6222_spike <- spikein_probe(GSE6222_raw)

GSE55092_deg <- RNAdegtable(GSE55092_rawdata)
GSE45436_deg <- RNAdegtable(GSE45436_rawdata)
GSE6222_deg <- RNAdegtable(GSE6222_raw)

deg_threshold <- max(GSE45436_deg[GSE45436_deg$actin3.actin5 < 3 & GSE45436_deg$gapdh3.gapdh5 < 1.25, "RNAdeg"])
little_deg <- rownames(GSE45436_deg)[GSE45436_deg$RNAdeg <= deg_threshold]

right_experiment <- rownames(GSE45436_spike)[GSE45436_spike$ployA_spikes == TRUE & GSE45436_spike$Bio_spikes == TRUE]


GSE40873_qc <- make_qc("GSE40873")
library(doParallel)








make_estimate(GSE6222_expression, "GSE6222")

library(sva)
test <- as.matrix(data.frame(GSE45436_final_expression, GSE6222_expression))
csif <- data.frame("sample" = colnames(test), "batch" = rep(c(1,2), c(40, 7)), "condition" = c(GSE45436_final_group, rep(c("Normal", "Tumor"), c(2, 5))))
modcombat <- model.matrix(~1, data = csif)
batch <- csif$batch
combat_test <- ComBat(dat = test, batch = batch, mod = modcombat, par.prior = T, prior.plots = T)

library(parallel)
cl.cores <- detectCores()
cl <- makeCluster(cl.cores)
system.time(GSE6222_spike <- clusterCall(cl, spikein_probe(GSE6222_rawdata)))
stopCluster(cl)