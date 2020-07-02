options(stringsAsFactors = F)

ssh        <- read.delim2("~/124.207.243.11/ssh_attack.txt", header = F, sep = ";", stringsAsFactors = F)
ssh_filter <- unique(ssh[,c(1,3)])
ssh_ip     <- unique(ssh_filter[,2])
ssh_ip     <- data.frame(ssh_ip, stringsAsFactors = F)


#delet space of ip character
cl     <- makeCluster(4)
ssh_ip <- parApply(cl, ssh_ip, 1, str_trim)
stopCluster(cl)
ssh_ip <- data.frame(ssh_ip[-605], stringsAsFactors = F)


cl <- makeForkCluster(6)
ssh$V3 <- parApply(cl, data.frame(ssh[, 3]), 1, str_trim)
stopCluster(cl)


library(RCurl)  #调用getURL()函数 
library(RJSONIO)  #调用fromJSON()函数 
library(stringr)
library(parallel)

Sinaurl <- function(ip) {
  paste("http://int.dpool.sina.com.cn/iplookup/iplookup.php?format=js&ip=", ip, 
    sep = "")
}  #sinaIP API 



fanxi <- function(aaa) {
  AA  <- NA
  BB  <- NA
  url <- NA
  cou <- NA
  pro <- NA
  cit <- NA
  
  ip <- NA  #定义初始值为0 
  for (i in 1:nrow(aaa)) {
    AA[i]  <- Sinaurl(aaa[i, 1])  #接口请求连接 
    url[i] <- getURL(AA[i])  #接口返回结果 
    BB[i]  <- strsplit(url[i], "=")
    BB[i]  <- gsub("^ ", "", BB[i][[1]][2])  #去掉首行空格 
    BB[i]  <- gsub(";", "", BB[i])  #去掉尾部分号 
    
    if(is.null(fromJSON(BB[[i]])[4:6]$country))
    {
      cou[i] <- "NULL"
      pro[i] <- "NULL"
      cit[i] <- "NULL"
      
    }
    else 
    {
      cou[i] <- fromJSON(BB[[i]])[4:6]$country  #提取国家 
      pro[i] <- fromJSON(BB[[i]])[4:6]$province  #提取省份 
      cit[i] <- fromJSON(BB[[i]])[4:6]$city  #提取城市 
    }
    ip[i] <- aaa[i, 1]
    #Sys.sleep(1)  #每次循环休眠1s 
  }
  return(data.frame(ip = ip, country = cou, province = pro, city = cit))  #汇总结果
}



# 定义结果输出列表
##MM <- list()
#n <- ceiling(nrow(Ip_yb)/100) - 1  #将原样本等分，除最后一份外，每份均含100个观测值 
n <- ceiling(nrow(ssh_ip)/100) - 1
#pb <- txtProgressBar(min = 0, max = n, style = 3)  #设置循环进度条 

library(foreach)
library(doParallel)
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
ssh_ip_translation <- foreach(i=1:n, .packages = c("RCurl", "RJSONIO"), .export = "fanxi", .combine = "rbind") %dopar% 
  {
  fanxi(data.frame(ssh_ip[(100i - 99):(100i), ], stringsAsFactors = F))  ##此处一定要注意添加stringsAsFactors=F，不然ip带不出来 
  }
#stopImplicitCluster() is can't work
stopCluster(cl)

colnames(ssh)[3] <- "ip"
ssh_result <- merge(data.frame(table(ssh$ip)), ssh_ip_translation, by = 1)
ssh_result <- unique(ssh_result)


library(ggplot2)
library(data.table)
library(ggtech)
test <- data.table(ssh_result)

country_count <- data.frame(table(ssh_result$country))
colnames(country_count)[1] <- "country"
test <- unique(merge(country_count,ssh_result, by="country"))
test <- data.table(test)

plot_data <- unique(test[, .(Freq.x, sum(Freq.y)), by=country])
plot_data <- plot_data[-1,]

ssh_plot  <- ggplot(plot_data, aes(country, Freq.x, colour = country)) + geom_count(aes(size= V2)) + theme_tech() + theme(axis.text.x = element_text(angle = 90)) +
  scale_size_area(max_size = 200) + ylab("IP") + guides(size = F, colour = F) + scale_y_continuous(trans = "log2")
 + theme(plot.subtitle = element_text(family = "sans"), 
    axis.text = element_text(size = 14), 
    axis.text.x = element_text(family = "mono"), 
    axis.text.y = element_text(family = "mono", 
        size = 25), plot.title = element_text(family = "sans"))
ggsave("./ssh_plot.jpeg", plot = ssh_plot, device = "jpeg", width = 40, height = 30, limitsize = F)
ggsave("./ssh_plot.jpeg", plot = ssh_plot, device = "jpeg", width = 15, height = 9)

