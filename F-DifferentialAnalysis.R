rm(list = ls())
dir.create("output")
##################################################
library('dplyr')
library(ggthemes)
library(ggplot2)
library(VennDiagram)
library(tidyverse)
library(VennDiagram)
library(limma)
library(RCircos)
library(pheatmap)
library(plyr)
library(officer)

#选择分析版本
#C-DifferentialAnalysis-Function-v2有报告，Function v1不生成报告
source("source/DifferentialAnalysis-Function-v2.R")

###功能说明：
#本程序先进行差异分析再与表型基因取交集并绘制火山图，热图，韦恩图和染色体定位图
#在验证集存在的情况下会将表型基因与验证集的基因取交集去掉不存在的表型基因

###输入参数说明：
#--dataset：可选“single”，“Two”。“single”模式只对单个数据集做差异分析，“Two”模式对两个数据集做差异分析取交集
#--exp：表达矩阵
#--exp_group：数据集分组
#--GSEID_ver：验证数据集前缀，无验证数据集填""
#--treat_name：疾病组的分组名
#--con_name：对照组的分组名
#--phenotype：表型基因路径
#--treat_col/con_col：疾病组/对照组颜色
#--heatmap_up_col/--heatmap_down_col：热图高表达/低表达颜色
#--volcano_width/volcano_height：火山图宽/高
#--heatmap_width/heatmap_height：热图宽/高

###input文件夹说明：
#2.1-差异分析.docx：用于替换的报告模板
#统一替换list.csv：程序需要读取该文件自动进行替换，需要修改第二列，顺序不可以改
#IRGs.csv：表型基因列表，可只有Gene Symbol列，第二列可为Genecard筛选分数
#Combined_Datasets_Group.csv：数据集分组
#Combined_Datasets_Matrix_norm.csv：表达矩阵
#GSE56815_Matrix_norm.csv：验证集表达矩阵

#单个数据集差异分析
hubgene<-Difference_analysis(dataset="single",exp="input/Combined_Datasets_Matrix_norm.csv",exp_group="input/Combined_Datasets_Group.csv",
                             GSEID_ver="",con_name="Control",treat_name="AP",treat_col="#A9907E",
                             con_col="#ABC4AA",phenotype="input/PRGs.csv",heatmap_up_col="#E64B35",heatmap_down_col="#4DBBD5",
                             volcano_width=8,volcano_height=7,heatmap_width=8,heatmap_height=8, chromosome_loc = F)

#两个数据集差异分析（暂时不用）
# hubgene<-Difference_analysis(dataset="Two",GSEID="Combined_Datasets",treat_name="OP",con_name="Control",GSEID2="GSE56815",
#                              GSEID_ver="",treat_col="#C99E8C",con_col="#465E65",
#                              phenotype="Ferroptosis",heatmap_up_col="#CD534C",heatmap_down_col="#7AA6DC",
#                              volcano_width=8,volcano_height=5,heatmap_width=8,heatmap_height=8)
