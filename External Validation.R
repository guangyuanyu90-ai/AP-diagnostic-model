#################################################ID转换、分组、数据标准化：
# Load necessary libraries
library(dplyr)
library(org.Hs.eg.db)  # For gene ID conversion

# Read the CSV file
data <- read.csv("GSE194331_HC_PAN_PANSEP_counts.csv")

# Convert GeneName from Ensembl to Symbol
data$GeneName <- mapIds(org.Hs.eg.db, keys=data$GeneName, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# Rename columns to remove hyphens and dots
colnames(data) <- gsub("[.-]", "", colnames(data))

# Save the modified data to a new CSV file
write.csv(data, "count.csv", row.names = FALSE)

# Load necessary libraries
library(dplyr)
library(DESeq2)

# Read the CSV file
data <- read.csv("count.csv")

# Convert GeneName from Ensembl to Symbol (if not already done)
# data$GeneName <- mapIds(org.Hs.eg.db, keys=data$GeneName, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# 处理重复基因名，按列（样本）求平均值
data_long <- data %>%
  pivot_longer(cols = -GeneName, names_to = "Sample", values_to = "Count") %>%
  group_by(GeneName, Sample) %>%
  summarise(Count = mean(Count)) %>%
  pivot_wider(names_from = "Sample", values_from = "Count")

# Standardization using DESeq2 method
# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = data_long[,-1], colData = data.frame(condition = colnames(data_long)[-1]), design = ~1)

# Perform DESeq2 normalization
dds <- estimateSizeFactors(dds)
data_normalized <- counts(dds, normalized = TRUE)

# Combine the normalized data with the GeneName
data_normalized <- cbind(GeneName = data_long$GeneName, data_normalized)

# Save the modified data to a new CSV file
write.csv(data_normalized, "/mnt/data/group.csv", row.names = FALSE)


################################################生成4个基因的数据集：
# Step 1: 读取数据
data_normalized <- read.csv("group_normalized—DESeq2.csv", check.names = FALSE)

# Step 2: 提取包含所需基因的行（假设这些基因的名称在 GeneName 列）
genes_of_interest <- c("VTN", "ANXA3", "RELA", "IQGAP1")

# 从 GeneName 列中选择包含这些基因的行
model_data <- data_normalized[data_normalized$GeneName %in% genes_of_interest, ]

# Step 3: 提取样本列，并将其与基因表达数据合并
# 将 GeneName 设置为行名，并去掉 GeneName 列
model_data <- model_data[, -1]  # 去除 GeneName 列
rownames(model_data) <- data_normalized$GeneName[data_normalized$GeneName %in% genes_of_interest]  # 设置行名为 GeneName

# Step 4: 创建分组信息：前 32 列为 Control 组，其余为 AP 组
group_data <- data.frame(
  V1 = colnames(data_normalized)[-1],  # 排除 GeneName 列
  Group = c(rep("Control", 32), rep("AP", ncol(data_normalized) - 33))  # 前32列为Control组，其余为AP组
)

# Step 5: 将分组信息添加到 model_data 中

colnames(model_data) <- group_data$Group  
# 确保样本顺序一致，添加分组信息
write.csv(model_data, "model_data.csv", row.names = TRUE)

##################################################先使用小鼠芯片的数据集进行AUC、Calibration诊断
# --- Step 1-7: 数据预处理 (保留您原有的转置和列名处理逻辑) ---
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(WGCNA)
library(tidyverse)
library(ggalluvial)
setwd("/media/desk16/iy5107/reviewer")


model_data <- read.csv("model_data.csv", check.names = FALSE)
model_data <- t(model_data)
model_data <- as.data.frame(model_data)
colnames(model_data) <- model_data[1, ]
model_data <- model_data[-1, ]

# --- 关键修正：结局变量 0/1 编码 ---
# 确保 Control = 0, AP = 1，避免 as.numeric 变成 1 和 2
genes <- c("ANXA3", "IQGAP1", "RELA", "VTN")
model_data[genes] <- lapply(model_data[genes], function(x) as.numeric(as.character(x)))
# --- Step 8: 应用文稿原始公式 (取代 glm 重新拟合) ---
model_data$Outcome_Num <- ifelse(grepl("AP", rownames(model_data), ignore.case = TRUE), 1, 0)
# 系数来源于 manuscript(4).docx Table 4
table(model_data$Outcome_Num)
intercept <- 2.062
coef_ANXA3 <- 0.21
coef_IQGAP1 <- 0.319
coef_RELA <- 2.754
coef_VTN <- -0.175


# --- 4. 计算预测概率 (不重新拟合模型) ---
# 计算线性预测值 LP
model_data$lp <- (model_data$ANXA3 * coef_ANXA3) + 
  (model_data$IQGAP1 * coef_IQGAP1) + 
  (model_data$RELA * coef_RELA) + 
  (model_data$VTN * coef_VTN) + 
  intercept

# 转换为预测概率 pred_probs
model_data$pred_probs <- 1 / (1 + exp(-model_data$lp))

library(pROC)
roc_obj <- roc(model_data$Outcome_Num, model_data$pred_probs, ci = TRUE)
cat("\n=== 严格外部验证结果 ===\n")
print(paste("AUC:", round(auc(roc_obj), 3)))

plot(roc_obj, legacy.axes = TRUE, col = "#D55E00", lwd = 3,
     main = "ROC Curve (Fixed Formula Validation)",
     print.auc = TRUE, print.auc.y = 0.4)

# --- 6. 生成校准图 (Calibration Plot) ---
library(ggplot2)
ggplot(model_data, aes(x = pred_probs, y = Outcome_Num)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  scale_x_continuous(limits = c(0, 1), name = "Predicted Probability") +
  scale_y_continuous(limits = c(0, 1), name = "Observed Proportion") +
  ggtitle("Calibration Plot (External Validation)") +
  theme_bw()

##############################################################
# ==============================================================================
# 1. 数据加载与预处理
# ==============================================================================
library(rms)
library(pROC)
library(dcurves)
library(dplyr)

# 读取数据
data <- read.csv("model_dataqq.csv")

# 清理列名 (去除可能存在的空格，确保变量名合法)
colnames(data) <- trimws(colnames(data)) # 去除首尾空格
# 如果第一列是无用的索引列(Unnamed)，删除它
if(colnames(data)[1] == "X" || colnames(data)[1] == "Unnamed..0") {
  data <- data[,-1]
}

# 确保 Group 是数值型 0 和 1
# 假设你的数据里已经是0/1。如果是字符，需要转换，例如:
# data$Group <- ifelse(data$Group == "Case", 1, 0)

# ==============================================================================
# 2. 环境设置 (rms 包特有要求)
# ==============================================================================
# 打包数据环境，这对 rms 包的函数运行至关重要
dd <- datadist(data)
options(datadist="dd")

# ==============================================================================
# 3. 构建模型 (Logistic Regression)
# ==============================================================================

# 模型 A: 全变量 (4个基因)
# lrm 是 rms 包的逻辑回归函数，功能比 glm 更强大
fit.full <- lrm(Group ~ ANXA3 + IQGAP1 + RELA + VTN, 
                data = data, x = TRUE, y = TRUE)

# 模型 B: 精简变量 (2个差异基因)
fit.reduced <- lrm(Group ~ ANXA3 + IQGAP1, 
                   data = data, x = TRUE, y = TRUE)

# 输出模型概况 (查看 AIC 和 P值)
cat("=== 4基因模型 AIC (越低越好) ===\n")
print(AIC(fit.full))
cat("\n=== 2基因模型 AIC (越低越好) ===\n")
print(AIC(fit.reduced)) 

# ==============================================================================
# 4. 绘制 ROC 曲线并对比
# ==============================================================================
# 计算 ROC 对象
roc.full <- roc(data$Group, predict(fit.full, type="fitted"))
roc.reduced <- roc(data$Group, predict(fit.reduced, type="fitted"))

# 绘图
plot(roc.full, col="red", main="ROC Curve Comparison")
plot(roc.reduced, add=TRUE, col="blue", lty=2)
legend("bottomright", 
       legend=c(paste0("4-Genes AUC=", round(auc(roc.full), 3)), 
                paste0("2-Genes AUC=", round(auc(roc.reduced), 3))),
       col=c("red", "blue"), lty=c(1, 2))

# 统计学检验：两个AUC是否有显著差异？
roc.test(roc.full, roc.reduced)

# ==============================================================================
# 5. 绘制校准曲线 (Calibration Curve) - 重点！
# ==============================================================================
# 使用 Bootstrap (B=1000) 进行内部验证，这是这种图的标准画法
# 这种方法画出来的线是平滑的，能看出真实的校准度
cal.full <- calibrate(fit.full, method="boot", B=1000)
cal.reduced <- calibrate(fit.reduced, method="boot", B=1000)

# 绘图
plot(cal.full, xlim=c(0,1), ylim=c(0,1), 
     xlab="Predicted Probability", ylab="Actual Probability",
     main="Calibration Curve: 4 Genes vs 2 Genes")
# 叠加 2基因模型的线
plot(cal.reduced, add=TRUE, col="blue")

legend("bottomright", 
       legend=c("4-Genes (Red)", "2-Genes (Blue)", "Ideal"), 
       col=c("black", "blue", "gray"), lty=c(1,1,2), pch=c(1,1,NA))

# 解读：如果蓝线(2基因)和黑线(4基因)重合，或者蓝线更靠近对角虚线，
# 说明2基因模型足够好，甚至更好。

# ==============================================================================
# 6. DCA 决策曲线 (Decision Curve Analysis)
# ==============================================================================
# dcurves 包绘制方法
# 转换数据格式以适应 dcurves
dca_data <- data %>%
  mutate(
    prob_4genes = predict(fit.full, type="fitted"),
    prob_2genes = predict(fit.reduced, type="fitted")
  )

dca_plot <- dca(Group ~ prob_4genes + prob_2genes, 
                data = dca_data,
                thresholds = seq(0, 0.9, by = 0.01),
                label = list(prob_4genes = "Model 4-Genes", prob_2genes = "Model 2-Genes"))

# 绘制 DCA
plot(dca_plot)

# ==============================================================================
# 7. (附加) 绘制列线图 (Nomogram) - 仅针对最终的 2基因模型
# ==============================================================================
# ==============================================================================
# 优化版列线图绘制代码
# ==============================================================================
# ==============================================================================
# R语言：生成高清无裁剪 Nomogram 的终极代码
# ==============================================================================
library(rms)

# 1. Ensure the environment and model are ready
dd <- datadist(data); options(datadist="dd")
fit.reduced <- lrm(Group ~ ANXA3 + IQGAP1, data = data, x = TRUE, y = TRUE)

# 2. Open PDF device for plotting (adjust width and height for better visibility)
pdf("My_Final_Nomogram.pdf", width = 12, height = 10)

# 3. Adjust margins to allow space for the "Risk of Disease" label
par(mar=c(14, 6, 4, 4))  # Increase bottom margin further to 14

# 4. Plot the nomogram with adjustments
nom <- nomogram(fit.reduced, 
                fun = plogis, 
                fun.at = c(0.05, 0.2, 0.5, 0.8, 0.95), 
                lp = FALSE, 
                funlabel = "Risk of Disease")

# Plot the nomogram
plot(nom, 
     xfrac = 0.3,       # Adjust the width of variable names
     cex.var = 1.2,     # Increase font size for variable names
     cex.axis = 1,      # Adjust axis labels size
     cex.funlabel = 1.5,  # Adjust font size of the label
     lwd = 2)           # Increase line width for better visibility

# 5. Close the PDF device to save the file
dev.off()

# The file "My_Final_Nomogram.pdf" will be saved in the current working directory

###################只绘制2个基因的ROC-DC
# ==============================================================================
# 最终版：2基因模型 (ANXA3 + IQGAP1) 综合评估图表
# 包含：ROC (带95% CI), Calibration, DCA
# ==============================================================================

library(rms)
library(pROC)
library(dcurves)
library(dplyr)
library(ggplot2) # DCA通常需要ggplot支持

# 1. 数据准备
data <- read.csv("model_dataqq.csv")
colnames(data) <- trimws(colnames(data))
if(colnames(data)[1] == "X" || colnames(data)[1] == "Unnamed..0") {
  data <- data[,-1]
}

# 2. RMS 环境设置 (必须步骤)
dd <- datadist(data)
options(datadist="dd")

# 3. 构建 2基因模型
# 使用 lrm (用于校准曲线)
fit.reduced <- lrm(Group ~ ANXA3 + IQGAP1, data = data, x = TRUE, y = TRUE)
# 计算预测概率 (用于ROC和DCA)
data$prob_2genes <- predict(fit.reduced, type="fitted")

# ==============================================================================
# 开始绘图：保存为高清 PDF
# ==============================================================================
pdf("Final_2Genes_Analysis.pdf", width = 10, height = 8)

# ------------------------------------------------------------------------------
# 图 1: ROC 曲线 (带 AUC 及 95% CI)
# ------------------------------------------------------------------------------
# 计算 ROC 和 CI
roc_obj <- roc(data$Group, data$prob_2genes)
ci_val <- ci(roc_obj) # 计算置信区间 [lower, median, upper]

# 格式化标签文字，例如: "AUC: 0.957 (95% CI: 0.910 - 0.992)"
auc_label <- sprintf("AUC: %.3f (95%% CI: %.3f - %.3f)", 
                     roc_obj$auc, ci_val[1], ci_val[3])

# 绘图
plot(roc_obj, 
     col = "#1c61b6",    # 蓝色
     lwd = 3,            # 线条加粗
     main = "ROC Curve: ANXA3 + IQGAP1",
     legacy.axes = TRUE, # 使用传统的 1-Specificity X轴
     print.auc = FALSE)  # 我们自己手动加标签，所以这里关掉默认的

# 添加美化的图例/文字
legend("bottomright", 
       legend = auc_label, 
       col = "#1c61b6", 
       lwd = 3, 
       bty = "n",       # 去掉图例边框
       cex = 1.2)       # 字体放大

# ------------------------------------------------------------------------------
# 图 2: 校准曲线 (Calibration Curve)
# ------------------------------------------------------------------------------
# 使用 Bootstrap 验证 (B=1000)
cal <- calibrate(fit.reduced, method="boot", B=1000)

plot(cal, 
     xlim = c(0,1), ylim = c(0,1),
     xlab = "Predicted Probability", 
     ylab = "Actual Probability",
     main = "Calibration Curve (Bootstrap Validation)",
     cex.axis = 1, cex.lab = 1.2)

# 添加图例解释
legend("bottomright", 
       legend = c("Ideal", "2-Genes Model", "Bias-corrected"), 
       lty = c(2, 1, 1), 
       col = c("gray", "black", "blue"),
       bty = "n")

# ------------------------------------------------------------------------------
# 图 3: DCA 决策曲线 (Decision Curve Analysis)
# ------------------------------------------------------------------------------
# 注意：dcurves 包通常使用 ggplot2 绘图
dca_res <- dca(Group ~ prob_2genes, 
               data = data,
               thresholds = seq(0, 0.9, by = 0.01),
               label = list(prob_2genes = "2-Genes Model"))

# 将 ggplot 对象打印到 PDF 中
p_dca <- plot(dca_res) + 
  ggtitle("Decision Curve Analysis") +
  theme_classic() +  # 使用经典简洁主题
  theme(legend.position = c(0.8, 0.8)) # 图例放右上角

print(p_dca)

# ==============================================================================
# 结束绘图
# ==============================================================================
dev.off() # 关闭设备，完成保存

cat("图表已生成完毕！请在您的文件夹中查找 'Final_2Genes_Analysis.pdf' 文件。")

################################################################################
# ==============================================================================
# 0. 环境准备与数据加载
# ==============================================================================
library(dplyr)
library(pROC)
library(caret)
library(ResourceSelection) # 用于 Hosmer-Lemeshow 检验
library(ggplot2)

# 读取数据
model_data <- read.csv("model_dataqq.csv")
colnames(model_data) <- trimws(colnames(model_data))
if(colnames(model_data)[1] == "X" || colnames(model_data)[1] == "Unnamed..0") {
  model_data <- model_data[,-1]
}

# 确保 Group 是数值型 0/1 (用于 glm 和 Brier 计算)
# 如果原始是字符 Case/Control，请先转换，这里假设已经是 0/1
# model_data$Group <- ifelse(model_data$Group == "Case", 1, 0) 

# ==============================================================================
# 1. 训练最终的 2基因模型 (ANXA3 + IQGAP1)
# ==============================================================================
# 核心修改：只使用筛选出的两个最佳基因
final_formula <- as.formula("Group ~ ANXA3 + IQGAP1")
model <- glm(final_formula, data = model_data, family = binomial())

# 获取线性预测值 (Linear Predictor) 和 预测概率
lp <- predict(model, type = "link")
pred_probs <- predict(model, type = "response")

# ==============================================================================
# 2. 生成 S1 - S5 表格
# ==============================================================================

# --- S1: 简易模型系数表 ---
model_coefficients <- coef(model)
s1_table <- data.frame(
  Gene = names(model_coefficients)[-1],  # 去掉截距
  Coefficient = model_coefficients[-1]
)
write.table(s1_table, "S1_Model_Coefficients.csv", sep = ",", quote = FALSE, row.names = FALSE)
cat("S1 表格已生成。\n")

# --- S2: 再校准参数 (Calibration Slope & Intercept) ---
# 计算校准截距 (Calibration-in-the-large): offset(lp)
cal_int_model <- glm(model_data$Group ~ offset(lp), family = binomial)
cal_intercept <- coef(cal_int_model)[1]

# 计算校准斜率 (Calibration Slope): y ~ lp
cal_slope_model <- glm(model_data$Group ~ lp, family = binomial)
cal_slope <- coef(cal_slope_model)[2]

# 计算 Brier Score
brier_score <- mean((pred_probs - model_data$Group)^2)

s2_table <- data.frame(
  Metric = c("Calibration_Intercept", "Calibration_Slope", "Brier_Score"),
  Value = c(round(cal_intercept, 4), round(cal_slope, 4), round(brier_score, 4))
)
write.table(s2_table, "S2_Recalibration.csv", sep = ",", quote = FALSE, row.names = FALSE)
cat("S2 表格已生成。\n")

# --- S3: AUC 总结 ---
roc_obj <- roc(model_data$Group, pred_probs, quiet = TRUE)
roc_ci <- ci(roc_obj)

s3_table <- data.frame(
  Model = "2-Genes (ANXA3 + IQGAP1)",
  AUC = round(roc_obj$auc, 3),
  CI_Lower = round(roc_ci[1], 3),
  CI_Upper = round(roc_ci[3], 3)
)
write.table(s3_table, "S3_AUC_Summary.csv", sep = ",", quote = FALSE, row.names = FALSE)
cat("S3 表格已生成。\n")

# --- S4: 详细模型统计 (含 OR 值) ---
model_summ <- summary(model)$coefficients
# 计算 OR 值 (Odds Ratio) 和 95% CI
or_vals <- exp(coef(model))
or_ci <- exp(confint(model))

s4_table <- data.frame(
  Term = rownames(model_summ),
  Coefficient = round(model_summ[, 1], 4),
  Std_Error = round(model_summ[, 2], 4),
  Z_Value = round(model_summ[, 3], 4),
  P_Value = format.pval(model_summ[, 4], eps = 0.001),
  OR = round(or_vals, 4),
  OR_Lower = round(or_ci[, 1], 4),
  OR_Upper = round(or_ci[, 2], 4)
)
write.table(s4_table, "S4_Model_Coefficients_Detailed.csv", sep = ",", quote = FALSE, row.names = FALSE)
cat("S4 表格已生成。\n")

# --- S5: 个体预测详情 ---
s5_table <- data.frame(
  Sample = rownames(model_data),
  Observed_Group = model_data$Group,
  Predicted_Prob = round(pred_probs, 4),
  Individual_Brier = round((model_data$Group - pred_probs)^2, 4)
)
write.table(s5_table, "S5_Model_Performance.csv", sep = ",", quote = FALSE, row.names = FALSE)
cat("S5 表格已生成。\n")

# ==============================================================================
# 3. 5x5 交叉验证 (Repeated CV) 与 校准图 (Calibration Plot)
# ==============================================================================

cat("\n正在进行 5x5 交叉验证...\n")

# 准备 Caret 数据：Caret 最好用因子型 (Factor) 作为 Outcome
cv_data <- model_data
cv_data$Group <- factor(cv_data$Group, levels = c(0, 1), labels = c("Control", "Case"))

ctrl <- trainControl(
  method = "repeatedcv", 
  number = 5,       # 5折
  repeats = 5,      # 重复5次
  classProbs = TRUE, 
  summaryFunction = twoClassSummary, 
  savePredictions = "final" # 保存预测值用于画图
)

# 训练 CV 模型 (仅 2基因)
set.seed(2025) # 设定种子保证可重复
fit_cv <- caret::train(
  Group ~ ANXA3 + IQGAP1, 
  data = cv_data, 
  method = "glm", 
  family = binomial(),
  trControl = ctrl, 
  metric = "ROC"
)

# 输出 CV 结果汇总
print(fit_cv)
write.table(fit_cv$results, "5x5_CV_Results_Summary.txt", sep = "\t", quote = FALSE)

# --- 绘制校准图 (Calibration Plot) ---

# 提取 CV 中的预测结果
cv_preds <- fit_cv$pred
# 将 Case/Control 转换回 0/1 以便计算均值
cv_preds$obs_num <- ifelse(cv_preds$obs == "Case", 1, 0)

# Hosmer-Lemeshow 分箱 (5 bins)
# 按照预测概率 (Case) 进行切分
cv_preds$bin <- cut(cv_preds$Case, 
                    breaks = quantile(cv_preds$Case, probs = seq(0, 1, by = 0.2)), 
                    include.lowest = TRUE, labels = FALSE)

# 计算每个 Bin 的平均预测概率和实际发生率
cal_plot_data <- cv_preds %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(Case),      # 预测概率均值
    mean_obs = mean(obs_num),    # 实际发生率
    n = n(),
    se = sqrt((mean_obs * (1 - mean_obs)) / n) # 标准误 (用于画误差棒)
  ) %>%
  ungroup()

# 计算 HL 检验 P 值 (用于显示在标题)
hl_res <- hoslem.test(cv_preds$obs_num, cv_preds$Case, g = 5)

# 使用 ggplot2 绘图
p_cal <- ggplot(cal_plot_data, aes(x = mean_pred, y = mean_obs)) +
  # 绘制对角参考线 (完美校准)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  # 绘制误差棒 (95% CI)
  geom_errorbar(aes(ymin = mean_obs - 1.96*se, ymax = mean_obs + 1.96*se), width = 0.02, color = "black") +
  # 绘制点和线
  geom_line(color = "#1c61b6", size = 1) +
  geom_point(size = 3, color = "#1c61b6") +
  # 坐标轴范围
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  # 标签和标题
  labs(
    title = paste0("Calibration Plot (5x5 CV, 2-Genes)\nHosmer-Lemeshow p = ", round(hl_res$p.value, 3)),
    subtitle = "Models: ANXA3 + IQGAP1",
    x = "Predicted Probability",
    y = "Observed Probability"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

# 保存图片
ggsave("Calibration_Plot_5x5CV.pdf", plot = p_cal, width = 6, height = 6)
print(p_cal)

cat("\n所有文件已生成：\n1. S1-S5 表格 (.csv)\n2. 5x5 CV 结果 (.txt)\n3. 校准曲线图 (Calibration_Plot_5x5CV.pdf)\n")

# ==============================================================================
# 绘制 Figure D: 分箱校准曲线 (Binned Calibration Plot)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(ResourceSelection) # 用于计算 HL test p-value

# 1. 准备数据 (确保使用最终的 2基因模型)
# -----------------------------------------------------------
data <- read.csv("model_dataqq.csv")
# 数据清理
if(colnames(data)[1] == "X" || colnames(data)[1] == "Unnamed..0") { data <- data[,-1] }
colnames(data) <- trimws(colnames(data))

# 训练模型
fit_final <- glm(Group ~ ANXA3 + IQGAP1, data = data, family = binomial())

# 获取预测概率
data$pred <- predict(fit_final, type = "response")

# 2. 数据分箱处理 (核心步骤)
# -----------------------------------------------------------
# 将预测概率分成 5 个区间 (Bins)，也可以改成 10 (g=10)
# cut() 函数把连续的概率切成几段
data$bin <- cut(data$pred, 
                breaks = quantile(data$pred, probs = seq(0, 1, 0.2)), # 5等分
                include.lowest = TRUE, 
                labels = FALSE)

# 计算每个箱子的：预测均值、实际发生率、标准误
cal_plot_data <- data %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),              # X轴: 平均预测概率
    mean_obs = mean(Group),              # Y轴: 实际发生率
    n = n(),                             # 样本量
    se = sqrt((mean_obs * (1 - mean_obs)) / n) # 标准误 (用于画误差棒)
  ) %>%
  ungroup()

# 3. 计算 Hosmer-Lemeshow 检验 (用于标题展示)
# -----------------------------------------------------------
hl_test <- hoslem.test(data$Group, data$pred, g = 5)
hl_p_value <- round(hl_test$p.value, 3)

# 4. 使用 ggplot2 绘图 (美化版)
# -----------------------------------------------------------
p <- ggplot(cal_plot_data, aes(x = mean_pred, y = mean_obs)) +
  # A. 绘制对角参考线 (完美预测线)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", size = 1) +
  
  # B. 绘制模型表现线 (连接各点)
  geom_line(color = "#1c61b6", size = 1) + 
  
  # C. 绘制误差棒 (95% CI) - 说明结果的稳健性
  geom_errorbar(aes(ymin = mean_obs - 1.96*se, ymax = mean_obs + 1.96*se), 
                width = 0.02, color = "black", alpha = 0.7) +
  
  # D. 绘制散点 (大小代表样本量)
  geom_point(aes(size = n), color = "#1c61b6", alpha = 0.8) +
  
  # E. 设置坐标轴范围 (0 到 1)
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  # F. 标题和标签
  labs(
    title = paste0("Calibration Plot (Hosmer-Lemeshow p = ", hl_p_value, ")"),
    subtitle = "Model: ANXA3 + IQGAP1",
    x = "Predicted Probability",
    y = "Observed Probability",
    size = "Sample Size" # 图例名称
  ) +
  
  # G. 主题美化 (经典学术风格)
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "bottom"
  )

# 5. 显示并保存
print(p)
ggsave("FigureD_Calibration_Binned.pdf", plot = p, width = 6, height = 6)


##############################################11个基因的箱型图：
# 加载必要的包
library(DESeq2)
library(ggplot2)
library(dplyr)

# 读取数据
count_data <- read.csv("group_normalized—DESeq2.csv")
head(count_data$GeneName, 20)  # 查看前20个基因名
count_data_clean <- count_data[!is.na(count_data$GeneName), ]
rownames(count_data_clean) <- count_data_clean$GeneName
count_data_clean <- count_data_clean[, -1]

group_info <- factor(c(rep("Control", 32), rep("AP", ncol(count_data_clean) - 32)))
genes_of_interest <- c("RELA", "IQGAP1", "ACTN4", "FLNA", "ANXA3", "VTN", "LCN2", "MPEG1", "PAH", "CA9", "PYHIN1")
vsd_genes <- count_data_clean[genes_of_interest, ]

# 将数据转换为长格式以便绘图
vsd_long <- as.data.frame(t(vsd_genes))
# 将 group_info 扩展为每个基因对应一个样本的组别
group_info <- c(rep("Control", 32), rep("AP", 87))

# 将group_info根据行数重复以匹配vsd_long的数据
vsd_long$Group <- c(rep("Control", 32), rep("AP", 87))

write.table(vsd_genes, "vsd_gene.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 1. 加载必要的绘图包

# 1. 加载包
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# 2. 读取数据
# 你的表格第一列似乎是样本名，各列是基因名 
df <- read.table("vsd_gene.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 3. 按照你的要求手动指定分组
# 前 32 个是 Control，后 87 个是 AP 
group_info <- c(rep("Control", 32), rep("AP", 87))
df$Group <- group_info

# 4. 进行 log2 转化 (使用 log2(x + 1) 防止 0 值报错)
# 我们只对前 11 列（基因列）进行计算 
df_log <- df
df_log<- log2(df + 1)

colnames(df_log) <- c(rep("Control", 32), rep("AP", 87))

write.table(df_log, "df_log.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# 1. 加载必要的库
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# 2. 读取数据 (假设文件名为 df_log.txt)
# 注意：你的文件第一行是分组名，后续行是基因名
lines <- readLines("df_log.txt")
groups <- unlist(strsplit(lines[1], "\t")) 
# 如果第一行只有 119 个元素，而基因行有 120 个，我们需要在 groups 前补一个空位
if(length(groups) == 119) {
  groups <- c("GeneName", groups)
}
# 2. 读取剩余的基因数据
# skip = 1 跳过第一行分组行
df_values <- read.table("df_log.txt", skip = 1, sep = "\t", header = FALSE, check.names = FALSE)

# 3. 命名并整理
colnames(df_values) <- groups
# 此时第一列已经是 GeneName，我们直接转换格式
df_long <- df_values %>%
  pivot_longer(cols = -GeneName, 
               names_to = "Group", 
               values_to = "Expression",
               names_repair = "minimal") %>%
  mutate(Expression = as.numeric(Expression)) # 确保数值类型正确

# 5. 绘图
ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.3) +
  facet_wrap(~ GeneName, scales = "free_y", ncol = 4) +
  scale_x_discrete(limits = c("Control", "AP")) +
  scale_fill_manual(values = c("Control" = "#4DBBD5FF", "AP" = "#E64B35FF")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")
) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold.italic", size = 10),
    legend.position = "none"
  ) +
  labs(x = "", y = "Log2 (Expression)")

# 1. 确保 Group 是因子且 Expression 是数值
df_long$Group <- factor(df_long$Group, levels = c("Control", "AP"))
df_long$Expression <- as.numeric(as.character(df_long$Expression))

# 2. 定义比较对象（解决 stat_compare_means 报错的关键）
my_comparisons <- list( c("Control", "AP") )

# 3. 绘图
ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.3) +
  facet_wrap(~ GeneName, scales = "free_y", ncol = 4) +
  # 关键修改：使用 label.y.npc 让 P 值位置自动适应每个子图
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format", 
                     label.x = 1.5,      # 位于 Control 和 AP 中间
                     label.y.npc = 0.9,  # 位于每个子图高度的 90% 处
                     size = 2.5) + 
  scale_fill_manual(values = c("Control" = "#4DBBD5FF", "AP" = "#E64B35FF")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold.italic", size = 10),
    legend.position = "none"
  ) +
  labs(x = "", y = "Log2 (Expression)")

###################四个基因判定严重度的ROC曲线绘制：
# 1. 加载必要的库
library(dplyr)
library(tidyr)

# 2. 读取 DESeq2 归一化后的数据集
# 假设你的文件在当前工作目录下
df_counts <- read.csv("group_normalized—DESeq2.csv", check.names = FALSE)



# 3. 定义你感兴趣的四个核心基因
target_genes <- c("ANXA3", "IQGAP1", "RELA", "VTN")


df_four_genes <- df_counts %>%
  filter(GeneName %in% target_genes)

# 5. 检查是否找齐了 4 个基因
if(nrow(df_four_genes) < length(target_genes)) {
  found_genes <- df_four_genes$GeneName
  missing_genes <- setdiff(target_genes, found_genes)
  warning("未能在 GeneName 列中找到以下基因: ", paste(missing_genes, collapse = ", "))
}


# 接你原来的代码：
df_four_genes <- df_counts %>%
  filter(GeneName %in% target_genes) %>%
  # 1. 将 GeneName 设为行名，这样它就不会参与 log2 计算
  column_to_rownames("GeneName") %>%
  # 2. 直接对整个数值矩阵进行 log2(x + 1) 转化
  # 此时数据框里全是数字，计算非常快
  mutate(across(everything(), ~log2(.x + 1))) %>%
  # 3. 转置：让样本变成“行”，基因变成“列”
  t() %>%
  as.data.frame()

# 此时 df_four_genes 的行名是 SampleID，列名是 4 个核心基因
# 预览一下结果
print(head(df_four_genes))

library(dplyr)
library(tibble)

# 1. 将行名转换成普通列 "SampleID"
df_analysis <- df_four_genes %>%
  rownames_to_column(var = "SampleID")

# 2. 此时数据类型如下：
# SampleID   ANXA3    IQGAP1    RELA      VTN
# HC025      7.461    13.093    8.869     0.000

library(dplyr)
library(tibble)
library(pROC)
library(ggplot2)

# 1. 读取分组信息
group_info <- read.csv("group.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(group_info) <- c("SampleID", "name", "RawGroup")

# 2. 重新定义严重程度分组（修正后的逻辑）
group_info <- group_info %>%
  mutate(Severity = case_when(
    # 匹配 Mild AP 或 MildAP
    grepl("Mild", RawGroup, ignore.case = TRUE) ~ "Mild AP",
    # 匹配 Moderately-severe 或 Severe AP
    grepl("severe", RawGroup, ignore.case = TRUE) ~ "Severe AP",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Severity)) # 确保只分析 AP 内部样本

# 3. 准备表达量数据并合并
df_analysis <- df_four_genes %>%
  rownames_to_column(var = "SampleID")

df_final <- inner_join(df_analysis, group_info, by = "SampleID")

# 4. 构建逻辑回归模型（联合 4 基因）
# Mild AP = 0, Severe AP (含中重度) = 1
df_final$Outcome <- ifelse(df_final$Severity == "Severe AP", 1, 0)

model <- glm(Outcome ~ ANXA3 + IQGAP1 + RELA + VTN, 
             data = df_final, family = binomial)

# 5. 计算预测概率
df_final$Prob <- predict(model, type = "response")

# 6. 计算 ROC 曲线
roc_obj <- roc(df_final$Outcome, df_final$Prob, ci=TRUE)

# 打印各组样本量，检查是否分类成功
print(table(df_final$Severity))
print(paste("AUC:", round(roc_obj$auc, 4)))

# 7. 绘制 ROC 曲线
ggroc(roc_obj, legacy.axes = TRUE, size = 1, color = "#21918c") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  annotate("text", x = 0.4, y = 0.2, 
           label = paste0("Joint 4-Gene AUC = ", round(roc_obj$auc, 3), 
                          "\n95% CI: ", round(roc_obj$ci[1], 3), "-", round(roc_obj$ci[3], 3)),
           size = 5) +
  theme_minimal() +
  labs(title = "Diagnostic Performance: Mild vs. (Mod-Sev & Severe) AP",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)")


#######################四个基因疾病严重度箱型图：
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 1. 读取并清洗分组信息
group_info <- read.csv("group.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(group_info) <- c("SampleID","name", "RawGroup")

# 使用方案 2：显式逻辑排除法进行分组
group_info <- group_info %>%
  mutate(Severity = case_when(
    # 匹配 Mild
    grepl("Mild", RawGroup, ignore.case = TRUE) ~ "Mild AP",
    
    # 包含 severe 但不包含 Moderately 的才是真正的重度
    grepl("Severe", RawGroup, ignore.case = TRUE) & 
      !grepl("Moderately", RawGroup, ignore.case = TRUE) ~ "SAP",
    
    # 包含 Moderately 的归为中重度
    grepl("Moderately", RawGroup, ignore.case = TRUE) ~ "MSAP",
    
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Severity))

# 设置因子水平顺序
group_info$Severity <- factor(group_info$Severity, levels = c("Mild AP", "MSAP", "SAP"))

# 2. 准备表达量数据 (假设 df_four_genes 已在内存中)
df_analysis <- df_four_genes %>%
  rownames_to_column(var = "SampleID")

# 3. 合并数据并转为长格式（方便 ggplot 分面绘图）
df_final <- inner_join(df_analysis, group_info, by = "SampleID")

df_long <- df_final %>%
  pivot_longer(cols = c(ANXA3, IQGAP1, RELA, VTN), 
               names_to = "Gene", 
               values_to = "Expression")

# 两两比较组合
my_comparisons <- list(
  c("Mild AP", "MSAP"),
  c("Mild AP", "SAP"),
  c("MSAP", "SAP")
)

p2 <- ggplot(df_long, aes(x = Severity, y = Expression, fill = Severity)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85) +
  geom_jitter(width = 0.18, size = 0.7, alpha = 0.35) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  # 先标注总体检验（可选）
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y.npc = 0.97,
    size = 3
  ) +
  # 再做两两比较（Wilcoxon = Mann-Whitney U）
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",         # 多重校正（推荐 BH/FDR）
    label = "p.signif",
    tip.length = 0.01,
    size = 3
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold.italic", size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(x = "", y = "Expression")

p2


######################严重度只分轻重的箱型图：
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 1. 重新定义二分类分组
group_info_binary <- group_info %>%
  mutate(Severity_Binary = case_when(
    grepl("Mild", RawGroup, ignore.case = TRUE) ~ "Mild AP",
    # 只要包含 severe (无论是 Moderately-severe 还是 Severe AP) 统统归为重度
    grepl("severe", RawGroup, ignore.case = TRUE) ~ "Severe AP",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Severity_Binary))

# 设置因子水平，固定左右顺序
group_info_binary$Severity_Binary <- factor(group_info_binary$Severity_Binary, 
                                            levels = c("Mild AP", "Severe AP"))

# 2. 合并表达量数据并转为长格式
df_analysis <- df_four_genes %>%
  rownames_to_column(var = "SampleID")

df_binary_final <- inner_join(df_analysis, group_info_binary, by = "SampleID")

df_long_binary <- df_binary_final %>%
  pivot_longer(cols = c(ANXA3, IQGAP1, RELA, VTN), 
               names_to = "Gene", 
               values_to = "Expression")

# 3. 绘制二分类箱线图
p_binary <- ggplot(df_long_binary, aes(x = Severity_Binary, y = Expression, fill = Severity_Binary)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4, color = "black") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  # 添加显著性检验 (P值或星号)
  stat_compare_means(method = "t.test", 
                     label = "p.signif", # 显示星号，若想显示具体数字请改为 "p.format"
                     label.x = 1.5,
                     label.y.npc = 0.9,
                     size = 5) +
  scale_fill_manual(values = c("Mild AP" = "#4DBBD5FF", "Severe AP" = "#E64B35FF")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold.italic", size = 12),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "Log2 (Expression)", 
       title = "Core Gene Expression: Mild vs. Severe/Mod-Sev AP")

print(p_binary)

# 4. 打印最终各组人数，确保合并正确
print("合并后的样本量统计：")
print(table(df_binary_final$Severity_Binary))

###################单独ANXA3的疾病严重度评分

library(dplyr)
library(tibble)
library(pROC)
library(ggplot2)

# 1. 读取分组信息 (确保只读取 2 列)
group_info <- read.csv("group.csv", header = FALSE, stringsAsFactors = FALSE)

# 检查列数并命名
if(ncol(group_info) == 2){
  colnames(group_info) <- c("SampleID", "RawGroup")
} else {
  colnames(group_info) <- c("SampleID", "name", "RawGroup") # 兼容你之前的写法
}

# 2. 重新定义严重程度分组
group_info <- group_info %>%
  mutate(Severity = case_when(
    grepl("Mild", RawGroup, ignore.case = TRUE) ~ "Mild AP",
    grepl("severe", RawGroup, ignore.case = TRUE) ~ "Severe AP",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Severity))

# 3. 准备表达量数据
df_analysis <- df_four_genes %>%
  rownames_to_column(var = "SampleID")

# 【关键点】合并数据并检查
df_final <- inner_join(df_analysis, group_info, by = "SampleID")

# --- 检查代码 ---
if (nrow(df_final) == 0) {
  stop("错误：合并后数据为空！请检查表达量数据的行名是否与 group.csv 的第一列匹配。")
}
if (length(unique(df_final$Severity)) < 2) {
  stop("错误：分组后只剩下一组数据（全是轻度或全是重度），无法计算 ROC。")
}
# ----------------

# 4. 定义结局变量
df_final$Outcome <- ifelse(df_final$Severity == "Severe AP", 1, 0)

# 5. 计算 ANXA3 的 ROC
# 加入 na.rm = TRUE 排除可能的空值
roc_anxa3 <- roc(df_final$Outcome, df_final$ANXA3, 
                 ci = TRUE, 
                 na.rm = TRUE)

# 6. 打印结果
print(table(df_final$Severity))
print(paste("ANXA3 Single Gene AUC:", round(roc_anxa3$auc, 4)))

# 7. 绘图
ggroc(roc_anxa3, legacy.axes = TRUE, size = 1, color = "#d95f02") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  annotate("text", x = 0.4, y = 0.2, 
           label = paste0("ANXA3 AUC = ", round(roc_anxa3$auc, 3)),
           size = 5) +
  theme_minimal() +
  labs(title = "Diagnostic Performance of ANXA3: Mild vs. Non-Mild AP")

#####################
# ==============================================================================
# code: Split-Sample Validation for Reviewer #3
# Purpose: To prove the model is not overfitted by training on 50% and testing on 50%
# ==============================================================================

# 1. 加载必要的包
if(!require(caret)) install.packages("caret")
if(!require(pROC))  install.packages("pROC")
library(caret)
library(pROC)

# 2. 读取数据 (请确保文件在工作目录下)
df <- read.csv("model_dataqq.csv", row.names = 1)

# 3. 数据整理
# 确保 Group 是因子变量 (0=Control, 1=AP)
df$Group <- as.factor(df$Group)
# 确保数值列是数值型 (有时候csv读取会变成字符)
df$ANXA3  <- as.numeric(df$ANXA3)
df$IQGAP1 <- as.numeric(df$IQGAP1)
df$RELA   <- as.numeric(df$RELA)
df$VTN    <- as.numeric(df$VTN)

# 检查一下数据概况
cat("总样本量:", nrow(df), "\n")
cat("分组情况:\n")
print(table(df$Group))

# ==============================================================================
# 核心步骤: 50% 训练 (Recalibration) / 50% 测试 (Validation)
# ==============================================================================

set.seed(123) # 设置种子，保证每次运行结果一致 (很重要!)

# 1. 创建拆分索引 (50% split)
trainIndex <- createDataPartition(df$Group, p = 0.5, list = FALSE, times = 1)

# 2. 生成 训练集 (Train) 和 测试集 (Test)
data_train <- df[trainIndex, ]
data_test  <- df[-trainIndex, ]

cat("\n[Split Summary]\n")
cat("Training Set (用于计算公式):", nrow(data_train), "samples\n")
cat("Testing Set  (用于验证AUC): ", nrow(data_test), "samples\n")

# 3. 在【训练集】上重新拟合模型 (Recalibration)
# 根据您之前的策略，人类模型主要依赖 ANXA3 和 IQGAP1 (2-gene panel)
# 如果您想用全部4个基因，把公式改为: Group ~ ANXA3 + IQGAP1 + RELA + VTN
fit_split <- glm(Group ~ ANXA3 + IQGAP1, data = data_train, family = binomial)

cat("\n[Model Recalibration]\n")
cat("Coefficients re-estimated on Training Set:\n")
print(coef(fit_split))

# 4. 在【测试集】上进行预测
# 注意：我们用训练集算出来的公式，去预测完全陌生的测试集
preds_test <- predict(fit_split, newdata = data_test, type = "response")

# 5. 计算 AUC (这是审稿人最关心的指标)
roc_obj <- roc(data_test$Group, preds_test, 
               levels = c("0", "1"), direction = "<", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ==============================================================================
# 结果输出
# ==============================================================================

cat("\n==========================================================\n")
cat("  FINAL RESULT FOR REVIEWER #3 (Split-Sample Validation)\n")
cat("==========================================================\n")
cat("  Test Set AUC =", round(auc_val, 3), "\n")
cat("==========================================================\n")

# 6. 绘图 (保存为PDF，留作证据)
pdf("Figure_Response_R3_SplitValidation.pdf", width = 6, height = 6)
plot(roc_obj, 
     col = "red", lwd = 2,
     main = paste0("Split-Sample Validation (50% Train / 50% Test)\nAUC = ", round(auc_val, 3)),
     legacy.axes = TRUE,
     xlab = "Specificity", ylab = "Sensitivity")
# 添加对角线
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

cat("图表已保存为: Figure_Response_R3_SplitValidation.pdf\n")

