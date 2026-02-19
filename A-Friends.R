rm(list = ls())
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("GOSemSim")
# install.packages("ggplot2")
library(org.Mm.eg.db)

gene <- data.table::fread("../7-SVM&LASSO/lasso_hubGenes.csv",data.table = F)
colnames(gene) <- "Gene"
gene <- clusterProfiler::bitr(gene$Gene,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
gene$ENTREZID <- as.character(gene$ENTREZID)
head(gene)
#用godata()函数来构建相应物种的Molecular Function本体的GO DATA
mf <- GOSemSim::godata('org.Mm.eg.db', ont="MF", computeIC = FALSE)
## preparing gene to GO mapping data...

#用godata()函数来构建相应物种的Cellular Component本体的GO DATA
cc <- GOSemSim::godata('org.Mm.eg.db', ont="CC", computeIC = FALSE)
## preparing gene to GO mapping data...

#用godata()函数来构建相应物种的Biological Process本体的GO DATA
bp <- GOSemSim::godata('org.Mm.eg.db', ont="BP", computeIC = FALSE)

#用mgeneSim来计算MF本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simmf <- GOSemSim::mgeneSim(gene$ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")

#用mgeneSim来计算CC本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simcc <- GOSemSim::mgeneSim(gene$ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")

#用mgeneSim来计算BP本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simbp <- GOSemSim::mgeneSim(gene$ENTREZID, semData = bp, measure = "Wang", drop = NULL, combine = "BMA")

#提取既有CC又有MF、BP的基因
comlist <- intersect(rownames(simcc), rownames(simmf))
comlist <- intersect(comlist, rownames(simbp))
simmf <- simmf[comlist, comlist]
simcc <- simcc[comlist, comlist]
simbp <- simbp[comlist, comlist]
row.names(gene) <- gene$ENTREZID
gene <- gene[comlist,]

#计算基因在MF本体和CC本体下的几何平均值，一个打分值同时包括基因的分子功能和细胞定位两个信息
fsim <- sqrt(simmf * simcc)
#或者计算基因在MF、CC、BP本体下的几何平均值
#fsim <- (simmf * simcc * simbp)^(1/3)

#将基因的名字由ENTREZID改为gene symbol，方便看懂。
colnames(fsim) = gene$SYMBOL
rownames(fsim) = gene$SYMBOL

#将基因自己和自己的相似度设为NA，方便接下来去掉。
for (i in 1:ncol(fsim)){
  fsim[i,i] <- NA
}

y <- reshape2::melt(fsim) #把宽格式数据转化成长格式，其实就是把正方形矩阵转成三列
y <- y[!is.na(y$value),] #删掉带NA的行

# 把每两个基因之间的相似度保存到文件，只需要保存第一列基因名和第三列数值
write.csv(y[,c(1,3)], "1-Friends_Output.csv", row.names = F)
y <- y[,c(1,3)]
head(y)
# Var1     value
# 2    APOB 0.5245798
# 3     LPL 0.4171678
# 4   PPARG 0.4404657
# 5     LEP 0.4764871
# 6    LEPR 0.4307784
# 7 SLC19A3 0.4901428
#计算每个基因跟其他基因相似度的平均值
y.mean <- aggregate(.~Var1,y,mean) 
m <- y.mean$value
names(m) <- y.mean$Var1
#按平均值给基因名排序，便于画图
y$Var1 <- factor(y$Var1, levels=names(sort(m)))

f <- function(y) {
  r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  r[3] <- mean(y)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(ggplot2)
p1 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  scale_fill_brewer(palette="Set2") + #配色
  guides(fill=FALSE) + #不显示图例
  
  stat_summary(fun.data= f, geom='boxplot') + 
  geom_hline(aes(yintercept=0.55), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细 
p1

# 保存到pdf文件
ggsave("1-Friends_Boxplot.pdf")
## Saving 4.34 x 4.99 in image

#画成云雨图
source("A-ViolinPlot-Function.R")
head(y)

# Var1     value
# 2    APOB 0.5245798
# 3     LPL 0.4171678
# 4   PPARG 0.4404657
# 5     LEP 0.4764871
# 6    LEPR 0.4307784
# 7 SLC19A3 0.4901428

unique(y$Var1)
## [1] "HSP90B1" "LSP1"    "ARHGAP4" "AHNAK"   "FLII"    "CWF19L1" "SIPA1"  
## [8] "FMNL1"

#计算每个分组的平均值
y.mean <- aggregate(.~Var1,y,mean) 
m <- y.mean$value
names(m) <- y.mean$Var1
#按平均值给分组排序，便于画图
y$Var1 <- factor(y$Var1, levels=names(sort(m))) 

p2 <- ggplot(y, aes(Var1, value, fill = Var1)) +
  scale_fill_brewer(palette="Set2") + #配色
  guides(fill=FALSE) +
  geom_flat_violin(position=position_nudge(x=.1)) +
  
  #分散不重叠的点图
  #geom_jitter(aes(color=Var1), width=.15) + guides(color=FALSE) +
  #堆叠的点图
  geom_dotplot(binaxis="y", stackdir="down", dotsize=.35) +
  stat_summary(fun.data= f, geom='boxplot') +
  # geom_boxplot(width=.1, position=position_nudge(x=.0)) +
  geom_hline(aes(yintercept=0.6), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细 

p2

ggsave("2-Friends_Raincloud2.pdf",width = 16, height = 5)

## Saving 4.34 x 4.99 in image
# #组个图，对比两种展示效果
# 
# library(cowplot)
# plot_grid(p1, p2, labels = c("A", "B"), align = "h")
# 
# ggsave("friends.pdf")
# 
# ## Saving 4.34 x 4.99 in image
# 
# sessionInfo()








