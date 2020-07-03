GSE45267_BC_plot <- ggplot(GSE45267_BC_data, aes(x=group, y=BC013423 , colour = group, shape = group))+
  geom_boxplot(width = 0.5)+scale_x_discrete(labels = c("Tumor", "Normal"), limits = c("T", "N"))+
  scale_y_continuous(breaks = c(2.2, 2.4, 2.6, 3), trans = "log2")+
  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
  theme_linedraw()+labs(title = "GSE45267-BC013423", x= "", y="Expression (log2)")+
  annotate("text",x=2.0,y=5.5,label = paste("Fold_change =", signif(GSE45267_dif["1568611_at", 1],4), sep = ""),hjust=0, size = 3)+
  annotate("text", x=2.0, y=5.2, label = paste("P.Val =", signif(GSE45267_dif["1568611_at", 4],4), sep = ""),hjust=0, size = 3)
ggsave(GSE45267_BC_plot, filename = "/home/yhy/liver_cancer/GSE45267_BC_plot.tiff", device = "tiff")



GSE55092_BC_plot <- ggplot(GSE55092_BC_data, aes(x=group, y= BC013423, colour = group, shape = group))+
  geom_boxplot(width = 0.5)+
  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
  scale_y_continuous(breaks = c(2.2, 2.4, 2.6, 3), trans = "log2")+scale_x_discrete(limits = c("Tumor", "Normal"))+
  #geom_violin(width = 1.1, trim = T)+
  theme_linedraw()+labs(title = "GSE55092-BC013423", x= "", y="Expression (log2)")+
  annotate("text",x=2.0,y=5.5,label = paste("Fold_change =", signif(GSE55092_dif["1568611_at", 1],4), sep = ""), hjust=0, size = 3)+
  annotate("text", x=2.0, y=5.2, label = paste("P.Val =", signif(GSE55092_dif["1568611_at", 4],4), sep = ""), hjust=0, size = 3)

ggsave(GSE55092_BC_plot, filename = "/home/yhy/liver_cancer/GSE55092_BC_plot.tiff", device = "tiff")


GSE45436_BC_plot <- ggplot(GSE45436_BC_data, aes(x=group, y= BC013423, colour = group, shape = group))+
  geom_boxplot(width = 0.5)+scale_x_discrete(labels = c("Tumor", "Normal"), limits = c("T", "N"))+
  scale_y_continuous(breaks = c(2.2, 2.4, 2.6, 3), trans = "log2")+
  geom_jitter(width = 0.3, alpha = 0.3, size = 2.5)+
  theme_linedraw()+labs(title = "GSE45436-BC013423", x= "", y="Expression (log2)")+
  annotate("text",x=2.0,y=5.5,label = paste("Fold_change =", signif(GSE45436_dif["1568611_at", 1],4), sep = ""),hjust=0, size = 3)+
  annotate("text", x=2.0, y=5.2, label = paste("P.Val =", signif(GSE45436_dif["1568611_at", 4],4), sep = ""),hjust=0, size = 3)
ggsave(filename = "/home/yhy/liver_cancer/GSE45436_BC_plot.tiff", plot = GSE45436_BC_plot, device = "tiff")


