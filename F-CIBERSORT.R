rm(list = ls())
dir.create("output")
library("Rmisc")
library("plyr")
library("Hmisc")
library("officer")
library("corrplot")
library('dplyr')

#选择分析版本
source("source/CIBERSORT-Function-V2.R")

###功能说明
#本程序先进行高低风险组的差异分析，随后绘制分组比较图，免疫细胞之间相关性热图，基因和免疫细胞的相关性气泡图

###输入参数说明：
#--exp：表达矩阵
#--exp_group：分组信息。
#treat_name/con_name：高风险组/低风险组名称，和分组信息中一致
#treat_col/con_col：高风险组/低风险组分组颜色
#heatmap_up_col/heatmap_down_col：热图上调/下调颜色
#immune_file：

###input文件夹说明：
#统一替换list.csv：程序需要读取该文件自动进行替换，需要修改第二列，顺序不可以改
#key_genes.csv：关键基因列表
#Combined_Disease_Matrix_norm.csv：表达矩阵
#LASSO_RiskGroup.csv：数据集分组

immu<-CIBERSORT_immu(exp="input/Combined_Disease_Matrix_norm.csv",exp_group = "input/LASSO_RiskGroup.txt",
                     immune_file = "source/mice.txt",
                     key_genes = 'input/lasso_hubGenes.csv',con_name="LowRisk",
                     treat_name="HighRisk",con_col="#8FBDD3",treat_col="#F97B22",
                     heatmap_up_col="#E64B35",heatmap_down_col="#4DBBD5")
