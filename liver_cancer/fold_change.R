library("limma")
library("lumi")

setwd("/home/yhy/liver_cancer/GSE81878/")
filenames_idat <- c("GSM2176944_NC1.idat", "GSM2176945_NC2.idat", "GSM2176946_NC3.idat", "GSM2176956_KDR1.idat", 
                    "GSM2176957_KDR2.idat", "GSM2176958_KDR3.idat")
#GSE81878_rawdata <- read.idat(idatfiles = filenames_idat, annotation = "Symbol")

idat2lumibatch <- function(filenames) {
      # filenames is a character vector of iDAT filenames
      require(illuminaio)
      require(lumi)
      idatlist = lapply(filenames,readIDAT)
      exprs = sapply(idatlist,function(x) {
          return(x$Quants$MeanBinData)})
      se.exprs = sapply(idatlist,function(x) {
          return(x$Quants$DevBinData/sqrt(x$Quants$NumGoodBeadsBinData))})
      beadNum = sapply(idatlist,function(x) {
          return(x$Quants$NumGoodBeadsBinData)})
#      detection = sapply(idatlist, function(x) {
#           return(x$Quants)})
      rownames(exprs)=rownames(se.exprs)=rownames(beadNum)=idatlist[[1]]$Quants$CodesBinData
      colnames(exprs)=colnames(se.exprs)=colnames(beadNum)=sapply(idatlist,function(x) {
          return(paste(x$Barcode,x$Section,sep="_"))})
      pd = data.frame(Sentrix=colnames(exprs))
      rownames(pd)=colnames(exprs)
      lb = new("LumiBatch",exprs=exprs,se.exprs=se.exprs,beadNum=beadNum,
                 phenoData=AnnotatedDataFrame(pd))
        
      return(lb)
  
}

#' Read idat files into a lumiBatch object
#'
#' Given a set of idat filenames and an annotation package,
#' returns a lumiBatch object with detection p-values calculated.
#'
#' @param filenames A character vector of the full path for each idat file to be read
#' @param annotation The name, as a character string (eg., "illuminaHumanv4.db"),
#' of the annotation package to be used for identifying bead types and for
#' inserting feature data.
#'
#' @export
#' 
idat2lumibatch_2 <- function(filenames,annotation) {
  require(illuminaio)
  require(lumi)
  #require(annotation, character.only = T)
  require(DBI)
  idatlist = lapply(filenames,readIDAT)
  exprs = sapply(idatlist,function(x) {
    return(x$Quants$MeanBinData)})
  se.exprs = sapply(idatlist,function(x) {
    return(x$Quants$DevBinData/sqrt(x$Quants$NumGoodBeadsBinData))})
  beadNum = sapply(idatlist,function(x) {
    return(x$Quants$NumGoodBeadsBinData)})
  rownames(exprs)=rownames(se.exprs)=rownames(beadNum)=idatlist[[1]]$Quants$CodesBinData
  colnames(exprs)=colnames(se.exprs)=colnames(beadNum)=sapply(idatlist,function(x) {
    return(paste(x$Barcode,x$Section,sep="_"))})
  pd = data.frame(Sentrix=colnames(exprs))
  rownames(pd)=colnames(exprs)
  conn = illuminaHumanv4_dbconn()
  tmp = dbGetQuery(conn,'select * from ExtraInfo')
  lb = new("LumiBatch",exprs=exprs,se.exprs=se.exprs,beadNum=beadNum,
           phenoData=AnnotatedDataFrame(pd))
  z = match(featureNames(lb),tmp$ArrayAddress)
  lb = lb[!is.na(z),]
  z = match(featureNames(lb),tmp$ArrayAddress)
  fData(lb)=tmp[z,]
  negprobes = which(fData(lb)$ReporterGroupName=='negative')
  #browser()
  assayDataElement(lb,'detection') = apply(exprs(lb),2,function(a) {
    return(illuminaPvalCalculation(a,negprobes))})
  return(lb)
}


#' Encapsulates a vectorized p-value calculation for Illumina microarrays
#'
#' The Illumina expression arrays use negative control probes for
#' determining the p-value of expression.  Refer to the Genome Studio
#' Gene Expression Module documentation for details.
#' 
#' @param values A vector of numeric values associated with a single sample
#' @param negProbeIdx An index defining the values that represent the negative
#' control probes
#'
#' @export
#' 
illuminaPvalCalculation = function(values,negProbeIdx) {
  ec = ecdf(values[negProbeIdx])
  return(1-ec(values))
}

#biocLite("illuminaHumanv4.db")
#biocLite("DBI")
library(illuminaHumanv4.db)
library(DBI)
GSE81878_rawdata <- idat2lumibatch_2(filenames_idat, illuminaHumanv4.db)
GSE81878_BG <- lumiB(GSE81878_rawdata, method = "bgAdjust")
GSE81878_VST <- lumiT(GSE81878_BG, method = "vst")
GSE81878_RSN <- lumiN(GSE81878_VST, method = "rsn")

GSE81878_QC <- lumiQ(GSE81878_RSN)
summary(GSE81878_QC, "QC")
plot(GSE81878_QC, what = "density")
plot(GSE81878_QC, what = "sampleRelation", method = "mds")

GSE81878_expression <- exprs(GSE81878_RSN)
presentCount <- detectionCall(GSE81878_rawdata)
GSE81878_expression_present <- GSE81878_expression[presentCount > 0, ]

array_probe <- toTable(illuminaHumanv4ARRAYADDRESS)
probe_symbol <- toTable(illuminaHumanv4SYMBOL)
colnames(array_probe) <- c("probe_id", "arrayaddress")
array_symbol <- merge(array_probe, probe_symbol, by.x = "probe_id")

#rownames(array_probe) <- array_probe$arrayaddress
GSE81878_expression_present <- as.data.frame(GSE81878_expression_present)
GSE81878_expression_present$arrayaddress <- rownames(GSE81878_expression_present)
GSE81878_expression_present <- merge(GSE81878_expression_present, array_symbol, by.x="arrayaddress")

#GSE81878_expression[rownames(array_probe), 7] <- array_probe[rownames(array_probe), 1]
#rownames(GSE81878_expression) <- GSE81878_expression[, 7]
#GSE81878_expression$V7 <- NULL

GSE81878_SYMBOL <- GSE81878_expression_present[, c(2,3,4,5,6,7,9)]
library("limma")
design <- model.matrix(~ -1 + factor(c(1,1,1,2,2,2)))
colnames(design) <- c("control", "RNAi")
fit <- lmFit(GSE81878_SYMBOL[, 1:6], design)
contrast.matrix <- makeContrasts(contrasts = "RNAi - control", levels = design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2, coef = "RNAi - control", n = nrow(fit2), lfc = log2(1.5))
dif <- dif[dif[, "P.Value"] < 0.01, ]
dif[rownames(dif), "SYMBOL"] <- GSE81878_SYMBOL[rownames(dif), 7] 

#fold_change <- read.csv("/home/yhy/liver_cancer/fold_change_report.csv")

SYMBOL_intersect <- intersect(fold_change$Gene.Symbol, dif$SYMBOL)
SYMBOL_intersect_detil <- dif[dif$SYMBOL %in% SYMBOL_intersect, ]
fold_change_intersect <- fold_change[fold_change$Gene.Symbol %in% SYMBOL_intersect, ]

UP_comm <- intersect(dif$SYMBOL[dif$logFC > 0], fold_change_intersect$Gene.Symbol[fold_change_intersect$log2.SK.762.SK.NC. > 0])
DOWN_comm <- intersect(dif$SYMBOL[dif$logFC < 0], fold_change_intersect$Gene.Symbol[fold_change_intersect$log2.SK.762.SK.NC. < 0])
UP_comm
DOWN_comm                                        

UP_dif <- intersect(dif$SYMBOL[dif$logFC > 0], fold_change_intersect$Gene.Symbol[fold_change_intersect$log2.SK.762.SK.NC. < 0])
DOWN_dif <- intersect(dif$SYMBOL[dif$logFC < 0], fold_change_intersect$Gene.Symbol[fold_change_intersect$log2.SK.762.SK.NC. > 0])
UP_dif
DOWN_dif


dif[dif$SYMBOL %in% "OLMALINC", ] # OLMALINC as down as hnRNPK

EntrezID_illumina <- toTable(illuminaHumanv4ENTREZID)
Entrez_array <- merge(array_probe, EntrezID_illumina, by = "probe_id")

arrayaddress_UP <- GSE81878_expression_present[rownames(dif)[dif$logFC > 0], "arrayaddress"]
arrayaddress_DOWN <- GSE81878_expression_present[rownames(dif)[dif$logFC < 0], "arrayaddress"]
arrayaddress_Total <- GSE81878_expression_present[rownames(dif), "arrayaddress"] 
illumina_Entrez <- unique(na.omit(Entrez_array$gene_id))


GO_analyse_illumina <- function(dif_arrayaddress){
  pearson_EG <- unique(Entrez_array$gene_id[Entrez_array$arrayaddress %in% dif_arrayaddress])
  params <- new("GOHyperGParams", geneIds = pearson_EG, universeGeneIds =illumina_Entrez, annotation = "illuminaHumanv4.db", ontology = "BP", pvalueCutoff = 0.001, conditional = F, testDirection = "over")
  hgOver <- hyperGTest(params)
  bp <- summary(hgOver)
  return(bp)
}

GSE81878_GO_UP <- GO_analyse_illumina(arrayaddress_UP)
GSE81878_GO_DOWN <- GO_analyse_illumina(arrayaddress_DOWN)
GSE81878_GO_Total <- GO_analyse_illumina(arrayaddress_Total)

for (i in c("GSE81878_GO_UP", "GSE81878_GO_DOWN", "GSE81878_GO_Total")) {
  write.csv(get(i), paste("/home/yhy/liver_cancer/", "Dif_", i, ".csv", sep = ""))
  
}







