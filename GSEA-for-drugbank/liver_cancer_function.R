##################################################################
#######
####### just get rawdata from celfile Dir.path, and change CDF
#######
###################################################################
get_rawdata1.0 <- function(GSE_name, not_qc = T, celfiles = NULL, custom_CDF = NULL){
          require("affy")
  
          path.cel <- paste0("~/liver_cancer/", GSE_name)
          if(is.null(celfiles)) celfiles <- list.celfiles(path.cel, full.names = T)
          else print("Used given celfiles \n")
          
          celfiles <- list.celfiles(path.cel, full.names = T)
          if(not_qc) print("Not qc \n")
          else {
            good_samples <- get(paste(GSE_name, "qc_result", sep = "_"))$good_sample
            number <- unlist(lapply(good_samples, FUN = function(x) grep(x, celfiles)))
            celfiles <- celfiles[number]
          }
          rawdata <- read.affybatch(celfiles)
          if(!is.null(custom_CDF)) {
            #All_custom_CDF <- c("ensg", "enst", "entrezg", "lncrna", "refseq")
            rawdata@cdfName <- paste(tolower(rawdata@cdfName), "hs", custom_CDF, "cdf", sep = "")
          }
          return(rawdata)
}

###################################################################
###
###              QC
###################################################################

make_qc_plm <- function(GSE_name){
              
  
              file_name <- paste("~/liver_cancer/", GSE_name, "/", sep = "")
              cel_file <- list.celfiles(file_name)
              setwd(file_name)
              rawdata <- read.affybatch(cel_file)
              qc_data <- qc(rawdata)
              qc_probe <- as.data.frame(qc_data@qc.probes)
              degradate <- rownames(qc_probe)[qc_probe$`AFFX-HSAC07/X00351_3_at`-qc_probe$`AFFX-HSAC07/X00351_5_at` > 3]
              PLM_data <- fitPLM(rawdata)
              NUSE_data <- as.data.frame(NUSE(PLM_data, type = "stats"))
              error <- colnames(NUSE_data)[NUSE_data[1,] > 1.02 | NUSE_data[1,] < 0.98]
              cel_file <- Reduce(setdiff, list(cel_file, degradate, error))
              return(list(cel_file = cel_file, degradate = degradate, error = error, qc = qc_data, PLM = PLM_data))
}

#####################################################################################
#
# RNA degradation result : make 3'/5' ratios of control probes and all probes together
#
#####################################################################################
RNAdegtable <- function(rawdata){
        require(simpleaffy)
        require(CLL)
        require(foreach)
        require(doParallel)
  
        cl <- makeCluster(2, type = "FORK")
        registerDoParallel(cl)
        result <- foreach(functions = list(qc, AffyRNAdeg)) %dopar% {
        do.call(functions, list(rawdata))
        }
        stopImplicitCluster()
        
        data_qc <- result[[1]]
        control_degratios <- as.data.frame(ratios(data_qc)[,c(1,3)]) #To get 3'/5' ratios of control probe: actin and GAPDH 
        RNAdeg <- result[[2]]
        RNAdeg <- as.data.frame(t(summaryAffyRNAdeg(RNAdeg)))[,1] # To calculate 3'/5' ratios of all probes  
        RNAdegtable <- data.frame(control_degratios, RNAdeg) # make them together to get RNAdegtable
        return(RNAdegtable)
}


###############################################################################
#
#make compare and PMA result of the two type of spike-in control probes together!
#
################################################################################

spikein_probe <- function(rawData){
        require("yaqcaffy")
        require("simpleaffy")
        require(parallel)
        require(foreach)
        require(doParallel)
  
        cl <- makeCluster(3, type = "FORK")
        registerDoParallel(cl)
        result <- foreach(functions = list(yaqc, detection.p.val, qc), .packages = c("yaqcaffy", "simpleaffy")) %dopar% {
          do.call(functions, list(rawData))
        }
        stopImplicitCluster()
        
        yack <- result[[1]]
        spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", #poly-A controls and only get 3' !
        rownames(yack@morespikes), ignore.case = TRUE),])
        sprep<-t(yack@morespikes[spnames,])
        # test for Lys < Phe < Thr < Dap
        t1 <- (sprep[,"lys3"] < sprep[,"phe3"] & sprep[,"phe3"] < sprep[,"thr3"] 
            & sprep[,"thr3"] < sprep[,"dap3"])       
        
        calls <- result[[2]]$call
        lys <- calls[rownames(calls)[grep("lys.*3",rownames(calls), ignore.case = TRUE)],]
        lys <- t(lys)  # get PMA of Lys, to make sure the lys probe that the minimum poly-A control probe present
        
        quality <- result[[3]]
        # test for BioB < BioC < BioD < CreX
        t2 <- (quality@spikes[,1]<quality@spikes[,2] 
                   & quality@spikes[,2]<quality@spikes[,3] 
                   & quality@spikes[,3]<quality@spikes[,4]) 
        BioB <- quality@bioBCalls # get PMA of BioB probe, to make sure the BioB probe that the minimum spike-in probe present
        spikes_result <- data.frame(list(ployA_spikes = t1, lys = lys, Bio_spikes = t2, BioB = BioB), stringsAsFactors = F)
        return(spikes_result)
}      


######################################################################
##
##Synthesize spike-in probes and RNA degradation to selest good sample
##
#######################################################################

good_sample <- function(spikes, RNAdeg){
          deg_threshold <- max(RNAdeg[RNAdeg$actin3.actin5 < 3 & RNAdeg$gapdh3.gapdh5 < 1.25, "RNAdeg"])# the max RNA degradation ratios of control probe
          little_deg <- rownames(RNAdeg)[RNAdeg$RNAdeg <= deg_threshold] 
          ## test for Lys < Phe < Thr < Dap and BioB < BioC < BioD < CreX, in addition probe of Lys and BioB have to present
          right_experiment <- rownames(spikes)[spikes$ployA_spikes == TRUE & spikes$Bio_spikes == TRUE & spikes$BioB == "P" & spikes$lys.AFFX.LysX.3_at == "P" & spikes$lys.AFFX.r2.Bs.lys.3_at == "P"]
          good_sample <- intersect(little_deg, right_experiment)
          return(list(little_deg = little_deg, right_experiment = right_experiment, good_sample = good_sample))
}

################################################################################################################
##
##Synthesize method of RNAdegtable and spikein_probe, and parallel computer each result of RNA degradation and spike-in probes
##Use function of good_sample to select good sample
##
#################################################################################################################

make_qc2.0 <- function(rawdata){
              require("yaqcaffy")
              require(simpleaffy)
              require(CLL)
              require(foreach)
              require(doParallel)
              ######################
              ## parallel compute ##
              ######################
              cl <- makeCluster(4, type = "FORK")
                registerDoParallel(cl)
                result <- foreach(functions = list(yaqc, detection.p.val, qc, AffyRNAdeg)) %dopar% {
                  do.call(functions, list(rawdata))
                }
              stopImplicitCluster()
              ##########################
              ###check spikes control###
              ##########################
              yack <- result[[1]]
              spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' !
                                                     rownames(yack@morespikes), ignore.case = TRUE),])
              sprep<-t(yack@morespikes[spnames,])
              # test for Lys < Phe < Thr < Dap
              t1 <- (sprep[,"lys3"] < sprep[,"phe3"] & sprep[,"phe3"] < sprep[,"thr3"] 
                     & sprep[,"thr3"] < sprep[,"dap3"])       
              ###################################
              #### check quality probe control###
              ###################################
              calls <- result[[2]]$call
              lys <- calls[rownames(calls)[grep("lys.*3",rownames(calls), ignore.case = TRUE)],]
              lys <- t(lys) 
              
              quality <- result[[3]]
              control_degratios <- as.data.frame(ratios(quality)[,c(1,3)])
              t2 <- (quality@spikes[,1]<quality@spikes[,2] 
                     & quality@spikes[,2]<quality@spikes[,3] 
                     & quality@spikes[,3]<quality@spikes[,4]) 
              BioB <- quality@bioBCalls
              ####################################
              ### check mean of rna degration ####
              ####################################
              RNAdeg <- result[[4]]
              RNAdeg <- as.data.frame(t(summaryAffyRNAdeg(RNAdeg)))[,1]
              RNAdeg <- data.frame(control_degratios, RNAdeg)
              spikes_result <- data.frame(list(ployA_spikes = t1, lys = lys, Bio_spikes = t2, BioB = BioB), stringsAsFactors = F)
              qc_result <- good_sample(spikes_result, RNAdeg)
              return(qc_result)
}

#################################################################################################
#
#Before do that fuction, get PMAtable first: PMAtable <- detection.p.val(rawdata)
#get each the threshold from given probes, decide the present threshold by this result
#
#################################################################################################

get_present_threshold <- function(PMAtable, probe_given){
                  PMAtable <- PMAtable$call
                  if(is.data.frame(probe_given)){
                  threshold <- apply(PMAtable[probe_given$probe,], 1, FUN = function(x) as.numeric(table(x, exclude = c("A", "M"))))
                  }
                  else{
                  threshold <- as.numeric(table(PMAtable[probe_given,], exclude = c("A", "M")))
                  }
                  threshold <- data.frame(threshold, probe_given)
                  return(threshold)
}

#############################################################
#
#filter all probes by the present threshold 
#
#############################################################

filter_probe <- function(PMAtable, present_threshold){
                PMAtable <- PMAtable$call
                present_probe <- rowSums(PMAtable == "P") >= present_threshold
                return(names(present_probe)[present_probe])
}

###########################################################################################################
#
# include the filter_probe and get_rawdata1.0 function, so I can get expression by given GSE_name. Before I 
# need qc result, make sure I make qc and have got GSE_qc_result
#
###########################################################################################################

get_expression1.0 <- function(GSE_name, method, present_threshold, PMAtable, custom_CDF = NULL, celfiles = NULL, not_qc = T){
          require("affy")
          require("gcrma")
          
          rawdata <- get_rawdata1.0(GSE_name, not_qc, celfiles, custom_CDF)
          switch(EXPR = method, "mas5" = {normalizedata <- mas5(rawdata)}, 
                 "rma" = {normalizedata <- rma(rawdata)}, 
                 "gcrma" = {normalizedata <- gcrma(rawdata)})
          expression <- as.data.frame(exprs(normalizedata))
          
          present_probes <- filter_probe(PMAtable, present_threshold)
          expression <- expression[present_probes,]
          if(!is.null(custom_CDF)) expression <- remove_at(expression)
          return(expression)
}

#########################################################################################
#
# This function include get_expression1.0 and parallel the three methods: rma, gcrma, mas5.
#
##########################################################################################

get_expression2.0 <- function(GSE_name, present_threshold, PMAtable, custom_CDF = NULL, celfiles = NULL, not_qc = T, custom_CDF_package = NULL){
      require("foreach")
      require("doParallel")
      cl <- makeCluster(3, type = "FORK")
      registerDoParallel(cl)
      expression <- foreach(method = list("rma", "gcrma", "mas5"), .export = c("get_expression1.0", "get_rawdata1.0", "filter_probe"),
                            .packages = custom_CDF_package) %dopar%{
        get_expression1.0(GSE_name, method, present_threshold, PMAtable, custom_CDF, celfiles, not_qc)
      }
      stopImplicitCluster()
      return(expression)
      
}
##############################################################################
##
## To get group or batch variable, factor_row is the group or batch lines of GEO_matrix
##
##############################################################################

get_group <- function(GEO_matrix, expression, factor_row, factor_name_list = NULL){ 
          order_number <- lapply(colnames(GEO_matrix), FUN = function(x){grep(x, colnames(expression))})
          group <- c()
          for(i in 1:length(order_number)){
            number <- order_number[[i]]
            if(length(number) == 0) next()
            else group[number] <- GEO_matrix[factor_row, i]
          }
          if(!is.null(factor_name_list)){
            group <- sub(factor_name_list[1], "Tumor", group)
            group <- sub(factor_name_list[2], "Normal", group)
          }
          return(group)
}

#############################################################################
##
##estimate the tumor pure, get estimate plot and return estimate score. results is in Document <- paste("~/liver_cancer/estimate_result/", GSE_Number, "/", sep = "")
##
#############################################################################

make_estimate <- function(expression, GSE_Number, database = "hgu133plus2.db"){
            require("estimate")
            require("annotate")
            require(database, character.only = T)
            
            expression <- get_unique_expression(expression, database = "hgu133plus2.db")
            Document <- paste("~/liver_cancer/estimate_result/", GSE_Number, "/", sep = "")
            filename <- paste(GSE_Number, "estimate", sep = "")
            dir.create(Document)
            write.table(expression, paste(Document, filename, ".txt", sep = ""), sep = "\t", quote = F)
            filterCommonGenes(input.f = paste(Document, filename, ".txt", sep = ""), output.f = paste(Document, filename, ".gct", sep = ""), id = "GeneSymbol")
            estimate_score_gct <- paste(Document, filename, "_score", ".gct", sep = "") 
            estimateScore(paste(Document, filename, ".gct", sep = ""), estimate_score_gct, platform="affymetrix")
            plotPurity(scores = estimate_score_gct, platform = "affymetrix", output.dir = paste(Document, "estimated_purity_plots", sep = ""))
            
            estimate_score <- read.delim(estimate_score_gct, row.names = 1, stringsAsFactors = F, skip = 2)
            estimate_score <- estimate_score[, -1]
            return(estimate_score)
}

###############################################################################
##
## get hierarchical clustering plot with color according to group and estimate of samples.
## make sure the order of group and estimate_score is same
##
###############################################################################

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
          
          jpeg(paste("~/liver_cancer/hc_plot/", name, ".jpeg", sep = ""), width = 1500, height = 1000)
          plot(hc_forplot)
          dev.off()
}

##########################################################################################
###
###remove the batch effect by ComBat. When expression is a list, I am merging many data sets
###
#############################################################################################

get_combat_expression <- function(expression, batch, group = NULL){
              require("sva")
              #what the fuck! if class(x)==data.frame, is.list(x)==Ture !!!
              if(!is.data.frame(expression)){
                  expression <- as.matrix(data.frame(unlist(expression)))
              }
              expression <- as.matrix(expression)
              
              if(!is.null(group)){
              csif <- data.frame("sample" = colnames(expression), "batch" = batch, "condition" = group)
              modcombat <- model.matrix(~1, data = csif)
              expression <- ComBat(dat = expression, batch = batch, mod = modcombat, par.prior = T, 
                                   prior.plots = T)
              }
              
              expression <- ComBat(dat = expression, batch = batch, par.prior = T, prior.plots = T)
              return(expression)
}

#######################################################################################################
#
#After filter probe expression and remove batch effect, I use MAX method to get unique expression and collpase to symbol or entrez.
#
#######################################################################################################

get_unique_expression <- function(expression, symbol = T, database = "hgu133plus2.db"){## database like "hgu133plus2.db"
                require(database, character.only = T)
                require("annotate")
                if(symbol){
                  col_name <- "SYMBOL"
                  probe_symbol <- as.data.frame(na.omit(getSYMBOL(rownames(expression), database)))
                  colnames(probe_symbol) <- col_name
                  expression <- data.frame(probe_symbol, expression[rownames(probe_symbol), ])
                  colnames(expression)[1] <- col_name
                }
                else{
                  col_name <- "ENTREZ"
                  probe_symbol <- as.data.frame(na.omit(getEG(rownames(expression), database)))
                  colnames(probe_symbol) <- col_name
                  expression <- data.frame(probe_symbol, expression[rownames(probe_symbol), ])
                  colnames(expression)[1] <- col_name
                }
                expression <- aggregate(as.formula(paste0(".~", col_name)), data = expression, max)
                rownames(expression) <- expression[, 1]
                expression <- expression[, -1]
                return(expression)
}

###########################################################################################
##
## make survival analyses, some thing wrong with my get_survival function, just use blow 'for' circulation
##
###########################################################################################

get_survival <- function(data, probe, gene, nmin = 10){
          require("survminer")
          require("survival")
          y <- Surv(LincRNA_survival$survival_day, LincRNA_survival$MO_status)
          probe_cutoff <- cutoff.survival(data, probe, nmin) 
          group <- LincRNA_survival[, probe] <= probe_cutoff$cutoff
          kmfit <- survfit(y~group)
          ggsurvplot(kmfit, pval = T, legend.title = paste("GSE40873-", gene, sep = ""), 
                     legend.labs = c("high_expression", "low_expression"), risk.table = T)
}

library(survminer)
library(survival)
library(progress)
#pb <- progress_bar$new(format = "survival ploting: (:spin) [:bar] :percent", total = nrow(probe_gene), clear = FALSE, width = 100)
#for (i in 1:nrow(probe_gene)) {
#      y <- Surv(LincRNA_survival$survival_day, LincRNA_survival$MO_status)
#     probe_cutoff <- cutoff.survival(LincRNA_survival, probe_gene$probe[i])
#     if (probe_cutoff$pvalue == 1) print(probe_gene$gene[i])
#     else 
#     group <- LincRNA_survival[, probe_gene$probe[i]] <= probe_cutoff$cutoff
#     kmfit <- survfit(y~group)
#     jpeg(paste("~/liver_cancer/GSE40873_survival/GSE40873", probe_gene$gene[i], Sys.Date(), "mas5", ".jpeg", sep = "-")) # 2017-01-03
#     print(ggsurvplot(kmfit, pval = T, legend.title = paste("GSE40873-", probe_gene$gene[i], sep = ""), legend.labs = c("high_expression", "low_expression"), risk.table = T))
#     dev.off()
#     pb$tick()
#                                                                                           

###################################################################################
##
## get fold change by Tumor/Normal
##
###################################################################################

get_flodchange <- function(group, eset, x = 1){
                require("limma")
                design2 <- model.matrix(~ -1 + factor(group))
                colnames(design2) <- c("Normal", "Tumor")
                fit <- lmFit(eset, design2)
                contrast.matrix <- makeContrasts(contrasts = "Tumor - Normal", levels = design2)
                fit1 <- contrasts.fit(fit, contrast.matrix)
                fit2 <- eBayes(fit1)
                dif <- topTable(fit2, coef = "Tumor - Normal", n = nrow(fit2), lfc = log2(x))
                #dif <- dif[dif[, "P.Value"] < 0.05, ]
                return(dif)
}

####################################################################
##
## make single probe differ express plot
##
####################################################################

get_plot <- function(probe_name, gene_name, expression, group, GSE_dif, GSE_number){
    require("ggplot2")
    GSE_data <- get_SYMBOL_data(probe_name, expression, group)
    GSE_plot <- ggplot(GSE_data, aes(y = SYMBOL, x = group , colour = group, shape = group))+
    geom_boxplot(width = 0.5)+scale_x_discrete(labels = c("Tumor", "Normal"), limits = c("Tumor", "Normal"))+
    geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
    theme_linedraw()+labs(title = paste(GSE_number, gene_name, sep = "-"), x= "", y="Expression (log2)")+
    annotate("text",x=0.5,y=5.5,label = paste("Fold_change =", signif(GSE_dif[probe_name, 1],4), sep = ""),hjust=0)+
    annotate("text", x=0.5, y=4.9, label = paste("adj.P.Val =", signif(GSE_dif[probe_name, 5],4), sep = ""),hjust=0)+
    annotate("text", x=0.5, y=5.2, label = paste("P.Val =", signif(GSE_dif[probe_name, 4],4), sep = ""),hjust=0)
    return(GSE_plot)
}

########################################################################
##
## make many probe differ express plot together, and save them to 
## paste("~/liver_cancer/dif_plot/", GSE_name, "/", sep = "")
##
########################################################################

get_plot2.0 <- function(expression, group, GSE_name, probe_gene){
      require("progress")
      dif <- get_flodchange(group, expression)
      path <- paste("~/liver_cancer/dif_plot/", GSE_name, "/", sep = "")
      if(!dir.exists(path)){dir.create(path)}
      
      pb <- progress_bar$new(format = paste("creating image -> ", path, "(:spin) [:bar] :percent", sep = ""),
        total = nrow(probe_gene), clear = FALSE, width = 100)
      for(i in 1:nrow(probe_gene)){
       differ_plot <- get_plot(probe_gene$probe[i], probe_gene$gene[i], expression, group, dif, GSE_name)
       ggsave(paste(path, probe_gene$gene[i], "_dif_", Sys.Date(),".jpeg", sep = ""), plot = differ_plot, device = "jpeg")     
       pb$tick()
      } 
}


################################################################################
##
## calculate the correlation of two probes
##
################################################################################

get_cor1.0 <- function(expression, probe1, probe2, group = NULL, met = "pearson"){
  if(is.null(group)){
      cor1 <- cor(as.numeric(expression[probe1, ]), as.numeric(expression[probe2, ]), method = met)
      cor1.p <- cor.test(as.numeric(expression[probe1, ]), as.numeric(expression[probe2, ]), method = met, alternative = "two.sided") 
      cor_result <- list(correlation = cor1, Pvalue = cor1.p)
  }
  else {
      group <- as.factor(group)
      level <- levels(group)
      cor_result <- list()
      for(i in 1:length(level)){
          level1 <- level[i]
          col1 <- grep(level1, group)
          print(level1)
          cor1 <- cor(as.numeric(expression[probe1, col1]), as.numeric(expression[probe2, col1]), method = met)
          cor1.p <- cor.test(as.numeric(expression[probe1, col1]), as.numeric(expression[probe2, col1]), method = met, alternative = "two.sided") 
          cor_result[[level1]] <-list(correlation = cor1, Pvalue = cor1.p) 
      }
  }
  return(cor_result)
}

#####################################################################################
##
## calculate the correlation of one probe with all probes of expression, and use p.adjust to adjust p value
##
#####################################################################################

get_cor2.0 <- function(expression, probe, group = NULL, met = "pearson", p.met ="bonferroni"){
        require("parallel")
        cl.cores <- detectCores()
        cl <- makeCluster(cl.cores, type = "FORK")
        if(is.null(group)){
            cor_result <- as.data.frame(parApply(cl, expression, 1, FUN = function(x){cor(as.numeric(expression[probe,]), as.numeric(x), method = met,use = "pair")}))
            colnames(cor_result) <- paste(probe, "cor", sep = "_")
            cor_result[, "pvalue"] <- parApply(cl, expression, 1, FUN = function(x){ y = cor.test(as.numeric(expression[probe,]), as.numeric(x), alternative = "two.sided", method = met);return(y$p.value)})
            cor_result$pvalue <- p.adjust(cor_result$pvalue, method = "bonferroni", nrow(cor_result))
        }
        else {
            group <- as.factor(group)
            level <- levels(group)
            cor_result <- list()
            
            for(i in 1:length(level)){
            level1 <- level[i]
            col1 <- grep(level1, group)
            print(level1)
            
            cor1 <- as.data.frame(parApply(cl, expression[, col1], 1, FUN = function(x){cor(as.numeric(expression[probe, col1]), as.numeric(x), method = met,use = "pair")}))
            colnames(cor1) <- paste(level1, probe, "cor", sep = "_")
            
            cor1[, "pvalue"] <- parApply(cl, expression[, col1], 1, FUN = function(x){ y = cor.test(as.numeric(expression[probe, col1]), as.numeric(x), alternative = "two.sided", method = met);return(y$p.value)})
            cor1$pvalue <- p.adjust(cor1$pvalue, method = "bonferroni", nrow(cor1))
            cor_result[[level1]] <- cor1
            }
        }
        stopCluster(cl)
        return(cor_result)
}

########################################################
####
#### GOstats, if there are group, select one of cor_result to make GO analyse. html report is saved in "~/liver_cancer/GOstats/"
####
########################################################

GO_analyse <- function(expression, cor_result, name, pos_threshold=0.3, neg_threshold=0.3, p=0.01, database = "hgu133plus2.db"){
            require("GOstats")
            require(database, character.only = T)
            gene_universe_expression <- get_unique_expression(expression, symbol = F)
            gene_universe <- rownames(gene_universe_expression)
            
            for(i in c("pos", "neg")){
              switch(EXPR = i, "pos" = {probe_select <- rownames(cor_result)[cor_result[,1] > pos_threshold & cor_result[, 2] < p]},
                               "neg" = {probe_select <- rownames(cor_result)[cor_result[,1] < neg_threshold & cor_result[, 2] < p]})
              expression_select <- expression[probe_select,]
              gene_select <- rownames(get_unique_expression(expression_select, symbol = F))
              params <- new("GOHyperGParams", geneIds = gene_select, 
                            universeGeneIds = gene_universe, annotation = database, ontology = "BP", 
                            pvalueCutoff = 0.001, conditional = F, testDirection = "over")
              hgOver <- hyperGTest(params)
              html_report <- paste(name, "GOstats", i, sep = "_")
              htmlReport(hgOver, file = paste0("~/liver_cancer/GOstats/", html_report, ".html"), label = html_report)
              
              bp <- paste("bp", i, sep = "_")
              assign(bp, summary(hgOver))
            }
            return(list(pos_GO = get("bp_pos"), neg_GO = get("bp_neg")))
}

####################
###
### GSEA functions!
###
#####################

############################################
##
## make the shell script that do GSEA.java saved in ~/GSEA-P-R/shell_script/, and make the PBS script to parallel computer saved in ~/GSEA-P-R/PBS_shell
### name is expression name, for example GSE****Tumor/Normal
###
############################################

make_script <- function(name, probe_gene, gmt_name = "c5.all.v5.2.symbols.gmt", GSEA_collapse = F){ 
          test <- readLines("~/GSEA-P-R/GSEA_test.sh")
          for (i in 1:nrow(probe_gene)) {
            script_name <- paste(name, "_", probe_gene$gene[i], sep = "")
            change_line <- sub("GSE45267_229090_at", script_name, test[2])
            change_line <- sub("GSE45267_Tumor", name, change_line)
            change_line <- sub("my_analysis_test", script_name, change_line)
            change_line <- sub("c5.all.v5.2.symbols.gmt", gmt_name, change_line) ##c5.bp.v5.2.symbols.gmt
            if(!GSEA_collapse){
            change_line <- sub("-chip DNAchip/HG_U133_Plus_2.chip ", "", change_line)
            change_line <- sub("-collapse true", "-collapse false", change_line)
            }
            filename <- paste("~/GSEA-P-R/shell_script/", script_name, ".sh", sep = "")
            write("#!/bin/bash", file = filename)
            write(change_line, file = filename, append = T)
            
            filename <- paste("~/GSEA-P-R/PBS_shell/PBS_", script_name, ".sh", sep = "")
            write("#!/bin/sh", file = filename)
            write("#PBS -N yhy", file = filename, append = T)
            write("#PBS -l mem=2gb,nodes=1:ppn=1", file = filename, append = T)
            write("#PBS -q 64G.q", file = filename, append = T)
            write("#PBS -V", file = filename, append = T)
            write("cd ~/yhy/GSEA-P-R/", file = filename, append = T)
            write(paste("sh shell_script/", script_name, ".sh", sep = ""), file = filename, append = T)
          }
}

############################################################
##
## get expression file of each single probe, saved in ~/GSEA-P-R/Datasets/
## make sure the probe present in the original_expression !!!
##
############################################################

phenotype_cls <- function(original_expression, name, probe_gene){
              for (i in 1:nrow(probe_gene)) {
                if(probe_gene$probe[i] %in% rownames(original_expression)){
                    filename <- paste("~/GSEA-P-R/Datasets/", name, "_", probe_gene$gene[i], ".cls", sep = "")
                    write("#numeric", file = filename, sep = "\t")
                    write(paste("#", probe_gene$gene[i], sep = ""), file = filename, sep = "\t", append = T)
                    write(unlist(original_expression[probe_gene$probe[i],1:ncol(original_expression)]), file = filename, ncolumns = ncol(original_expression), sep = "\t", append = T)
                }
                else cat(paste(probe_gene$probe[i], "is not present in expression"), "\n")
              }
}

phenotype_cls2 <- function(group_info, name){
                  filename <- paste0("~/GSEA-P-R/Datasets/", name, ".cls")
                  first <- ifelse(group_info[1]=="T", "Tumor", "Normal")
                  second <- ifelse(group_info[1]=="T", "Normal", "Tumor")
                  write(paste(length(group_info), "2", "1", sep = "\t"),file = filename, sep = "\t")
                  write(paste("#", first, second, sep = "\t"), file = filename, sep = "\t", append = T)
                  write(group_info, file = filename, sep = "\t", append = T, ncolumns = length(group_info))
}
##############################################################
##
## save the expression file in ~/GSEA-P-R/Datasets/
##
##############################################################

expres_gct <- function(expression, name){
          expression <- na.omit(expression)
          expression_GSEA <- cbind(Description = "na", expression)
          filename_exprs <- paste("~/GSEA-P-R/Datasets/", name, ".gct", sep = "")
          write("#1.2", file = filename_exprs)
          write(as.character(c(nrow(expression), ncol(expression))), file = filename_exprs, ncolumns = 2, sep = "\t", append = T)
          write(c("NAME", colnames(expression_GSEA)), file = filename_exprs, ncolumns = ncol(expression_GSEA)+1, sep = "\t", append = T)
          write.table(expression_GSEA, file = filename_exprs, sep = "\t", append = T, quote = F, col.names = F, row.names = T)
}

################################################################################################
#####
##### before do this function, make sure do sh clear_files first!
##### To make all files of GSEA need, expression_list is list(original expression, unique symbol expression, custom CDF or unique entrezg, .......)
#####
################################################################################################

make_all_GSEA_files <- function(expression_list, GSE_name, probe_gene, sample_group = NULL, gmt_name = "c5.all.v5.2.symbols.gmt", GSEA_collapse = F){
                require("doParallel")
                require("foreach")
                cl <- makeCluster(length(expression_list)-1, type = "FORK")
                registerDoParallel(cl)
                
                foreach(j = as.list(2:length(expression_list)), .export = c("make_script", "expres_gct", "phenotype_cls")) %dopar%{
                  if(!is.null(expression_list[[j]])){
                    expression <- expression_list[[j]]
                    expression_name <- c(NA, "symbol", "entrez")
                    expression_name <- expression_name[j]
                    if(j == 3) gmt_name <- "c5.all.v5.2.entrez.gmt"
                    if(!is.null(sample_group)){
                        level <- levels(as.factor(sample_group))
                        for(k in level){
                          level_expression <- expression[,grep(k, sample_group)]
                          level_name <- paste(GSE_name, expression_name, k, sep = "_")
                          expres_gct(level_expression, level_name)
                          make_script(level_name, probe_gene, gmt_name, GSEA_collapse = GSEA_collapse)
                          phenotype_cls(expression_list[[1]][,grep(k, sample_group)], level_name, probe_gene)
                        }
                    }
                    else {
                          GSE_name <- paste(GSE_name, expression_name, sep = "_")
                          expres_gct(expression, GSE_name)
                          make_script(GSE_name, probe_gene, gmt_name, GSEA_collapse = GSEA_collapse)
                          phenotype_cls(expression_list[[1]], GSE_name, probe_gene)
                    } 
                  }
                }
                  stopImplicitCluster()
}



















######################################################
#
# little function
#
######################################################
######################################################
#get the last n letter of x
######################################################

str_right <- function(x, n){
              substr(x, nchar(x)-n+1, nchar(x)-n+1)
}
########################################################################################
########################when used custom CDF, change *****_at to real name #############
########################################################################################

remove_at <- function(expression){
              control_rows <- grep("affx",tolower(rownames(expression)))
              expression <- expression[-control_rows,]
              rownames(expression) <- substring(rownames(expression), 1, nchar(rownames(expression))-3)
              return(expression)
}

#########################################################
#simply make hierarchical clustering
#########################################################

get_hc <- function(expression, method="com"){
  expression_cor <- cor(expression)
  expression_hc <- hclust(as.dist(1-expression_cor), method)
  return(expression_hc)
}

#######################################################
##cut off survival from Hu ShuoFeng, I change it
#######################################################

cutoff.survival <- function(data, probe, nmin=10) {
              require("survival")
              data <- data[order(data[,probe], decreasing = F),]
              marker <- data[,probe]
              time <- data$survival_day
              event <- data$MO_status
              n <- length(marker)
              nlow <- (nmin+1):(n-nmin)
              y <- Surv(time, event)
              p_min <- 1
              cutoff <- 0
              p_vector <- vector()
              for (i in nlow) {
                # j <- i-nmin+1   
                if (marker[i] != marker[i+1]) {
                  x <- c(rep(0, i), rep(1, n-i))
                  model <- summary(coxph(y ~ x))
                  # coef <- model$coefficients
                  p <- model[["sctest"]]["pvalue"]
                }else{p <- 2}
                p_vector <-c(p_vector,p)
                cutoff <- ifelse(p < p_min,i,cutoff)
                p_min <- ifelse(p < p_min,p,p_min)
              }
              p_vector <- c(rep(2,10),p_vector,rep(2,10))
              ss <- cbind(marker,p_vector)
              return(list(cutoff = marker[cutoff],pvalue = p_min))
              # return(marker[cutoff])
}

########################################################
## select probe expression and adjust data for make differ express plot
########################################################

get_SYMBOL_data <- function(probe_id , expression, data_group){
  GSE_data <- expression[probe_id,]
  GSE_data["group", ] <- data_group
  rownames(GSE_data)[1] <- "SYMBOL" #SYMBOL is character
  GSE_data <- t(GSE_data)
  GSE_data <- as.data.frame(GSE_data)
  GSE_data[["SYMBOL"]] <- as.numeric(as.character(GSE_data[["SYMBOL"]]))
  return(GSE_data)
}


####################################################################################
#
# How to install "Rmpi" for parallel, but I can't use it, I don't know why
#install.packages("Rmpi", configure.args = c("--with-mpi=/usr/local/lib", "--with_Rmpi-type=OPENMPI")) before that have to add "export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:}/usr/lib64/openmpi/lib/" to ~/.zshrc and resource it
#sudo sh -c 'echo /usr/local/lib > /etc/ld.so.conf.d/openmpi.conf; ldconfig'
#
######################################################################################
