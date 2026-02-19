rm(list = ls())
##################################################
library("limma")
library('dplyr')
###数据读取######
## ES差异分析
ES_GSE_mat <- as.data.frame(data.table::fread("../2-RawData/1-GEO/output/Merged_Dataset/Combined_Disease_Matrix_norm.csv"))
rownames(ES_GSE_mat) <- ES_GSE_mat$V1
ES_GSE_mat <- ES_GSE_mat[,-1]
ES_sampleinfo <- as.data.frame(data.table::fread("../7-SVM&LASSO/LASSO_RiskGroup.txt",header = F))
colnames(ES_sampleinfo) <- c('geo_accession', 'Group')
ES_GSE_mat <- ES_GSE_mat[,ES_sampleinfo$geo_accession]

##输出表型基因分组比较图
genes <- read.table("../7-SVM&LASSO/lasso_HubGenes.csv", header = T, sep = ',')
genelist <- ''
for(i in genes$x){
  genelist <- paste0(genelist, '，', i)
}
genelist

ES_boxpgse <- ES_GSE_mat[genes$x,]
all(ES_sampleinfo$geo_accession == colnames(ES_GSE_mat))

library(ggplot2)
library(corrgram)
library(ggthemes)
library(ggpubr)

# rownames(GRDEGs_matrix) <- NULL
ESRDEGs_long <- ES_boxpgse %>%
  # tibble::column_to_rownames("V1") %>%
  t() %>% #行列置换把Gene从列换成行
  data.frame() %>%
  dplyr::mutate(Sample = rownames(.)) %>% #把行名添加到Sample列，用于数据合并
  merge(ES_sampleinfo, by.x = "Sample", by.y = "geo_accession")  %>%
  tidyr::gather(key = geneid, value, -c(Sample, Group))
ESRDEGs_long$Group <- factor(ESRDEGs_long$Group, levels = c("LowRisk", "HighRisk"))

ES_barplot<-ggplot(ESRDEGs_long,aes(x=geneid,y=value,fill=Group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#8FBDD3","#F97B22"))+
  theme(panel.grid = element_blank())+
  stat_compare_means(method	= "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,color="black"))+
  theme(legend.position = 'top')+xlab('')+ylab('Gene Expression')+labs(fill='Group')+
  theme(text = element_text(size = 14))
ES_barplot
ggsave(ES_barplot,file="1-Boxplot_Risk.pdf",width = 16.6, height = 5)


##差异基因ROC
ES_heatmap_c1 <- t(ES_GSE_mat[genes$x,])
all(rownames(ES_heatmap_c1) == ES_sampleinfo$geo_accession)
# write.csv(ES_heatmap_c1, 'heatmap_input.csv', row.names = F)
ES_heatmap_c1 <- cbind(ES_heatmap_c1, ES_sampleinfo$Group)
colnames(ES_heatmap_c1)[5] <- 'Group'
for(i in 1:4){
  gname <- colnames(ES_heatmap_c1)[i]
  ms_roc <- ES_heatmap_c1[, c(5, i)]
  write.csv(ms_roc, paste0('ROC-', gname, '.csv'), row.names = F)
}
#write.csv(heatmap_c1, 'ROC.csv', row.names = F)
# ms <- tapply(1:8, ceiling((1:8)/3), function(i) ES_heatmap_c1[,i])
# for(i in 1:3){
#   ms_roc <- cbind(ES_sampleinfo$Group, ms[[i]])
#   colnames(ms_roc)[1] <- 'Group'
#   write.csv(ms_roc, paste0('ROC-', i, '.csv'), row.names = F)
# }
# 
# write.csv(ES_heatmap_c1, 'heatmap_input.csv', row.names = F)

