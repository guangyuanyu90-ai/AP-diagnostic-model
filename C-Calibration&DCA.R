#rm(list = ls())
library(ggplot2)
library(ggthemes)
library(ggpubr)

###整理输入文件######
data <- data.table::fread("../2-RawData/1-GEO/output/Merged_Dataset/Combined_Datasets_Matrix_norm.csv",data.table = F)
rownames(data) <- data$V1
data<-data[,-1]
data<-as.data.frame(t(data[,]))
group <- as.data.frame(data.table::fread("../2-RawData/1-GEO/output/Merged_Dataset/Combined_Datasets_Group.csv",header = T))
colnames(group) <- c("geo_accession", "Group")
table(group$Group)
data$cla<-c(rep('Control',17),rep('AP',53))
colnames(data)
phon<-read.table("../7-SVM&LASSO/LASSO_RiskScore.txt",header = T)
phon <- phon[group$geo_accession,]
phon$cla<-c(rep('Control',17),rep('AP',53))#再次确认，可以不运行
all(rownames(phon) == group$geo_accession)
library(rms)

###画LASSO模型的Nomogram图######
dc<-datadist(phon)
options(datadist="dc")
fit <- lrm(cla~RiskScore, data=phon, x=T, y=T,tol=1e-9,maxit=1000)
nom <- nomogram(fit,fun=plogis,fun.at =c(0.5,1) ,funlabel=c("Risk of Sarcopenia"))
#pdf(width = 8.1,height = 5,file = "Nomogram_1.pdf",onefile = F)
plot(nom,cex.axis=1.5,cex.var=1.5)
dev.off()

###画LASSO模型中5个基因的nomogram图######
genelist <- ""
for(i in colnames(phon)){
  genelist <- paste0(genelist, '+', i)
}
genelist
colnames(phon)
formula <- as.formula(paste('cla~', paste("`",colnames(phon)[1:(ncol(phon)-2)],"`",sep = "",collapse = "+")))
#phon1 <- colnames(phon)[1:20]
fit <- lrm(formula, 
           data=phon, x=T, y=T,tol=1e-9,maxit=1000)
dc<-datadist(phon)
options(datadist="dc")
nom <- nomogram(fit,fun=plogis,fun.at =c(0.5,1) ,funlabel=c("Risk of Sarcopenia"))
pdf(width = 8.1,height = 7,file = "2-Nomogram.pdf",onefile = F)
plot(nom,cex.axis=1.2,cex.var=1.5)
dev.off()

###画nomogram的校正曲线图######
#使用0代表分组为Normal的样本，使用1代表分组为Disease的样本
phon$cla<-c(rep('0',17),rep('1',53))
phon$cla <- phon$cla%>% as.numeric()
glm <- glm(formula,
           data=phon,family = binomial(link = "logit"),control=list(maxit=1000))
# glm <- glm(cla~RiskScore,
#            data=phon,family = binomial(link = "logit"),control=list(maxit=1000))
P1 <- predict(glm,type = 'response')
pdf(file = "3-CalibrationCurve.pdf",width = 9.01,height = 8.5)
val.prob(P1,phon$cla,cex = 1.5)
dev.off()

###画DCA图######
library(ggDCA)
library("rmda")
# simple<- decision_curve(cla~RiskScore,data = phon, family = binomial(link ='logit'),
                         # thresholds= seq(0,1, by = 0.01),
                         # confidence.intervals =0.95,study.design = 'case-control',
                         # population.prevalence = 0.3)
complex<-decision_curve(formula,data = phon,
               family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
               confidence.intervals= 0.95,study.design = 'case-control',
               population.prevalence= 0.3)
# List<- list(simple,complex)

#画lasso模型和lasso模型20个基因的DCA图
#pdf(width = 12,height = 8,file = "DCA_All.pdf")
# plot_decision_curve(List,curve.names= c('simple','complex'),
#                     cost.benefit.axis =FALSE,col = c("#83a78d", "#f0a274"),
#                     confidence.intervals =FALSE,standardize = FALSE)
# dev.off()

#单独画LASSO回归DCA图
#pdf(width = 12,height = 7,file = "6-DCA.pdf")
# plot_decision_curve(simple,curve.names= c('RiskScore'),
#                     cost.benefit.axis =FALSE,col = c("#f0a274"),
#                     confidence.intervals =FALSE,standardize = FALSE,xlim=c(0,0.6))
# dev.off()

#单独画3个gene的LASSO回归DCA图
pdf(width = 8.1,height = 8,file = "4-DCA.pdf")
plot_decision_curve(complex,curve.names= c('Model'),
                    cost.benefit.axis =FALSE,col = c("#eea2a4"),
                    confidence.intervals =FALSE,standardize = FALSE,xlim=c(0,0.75))
dev.off()
###END######
