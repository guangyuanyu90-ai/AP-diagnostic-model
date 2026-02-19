## =========================================================
## Internal validation (4-gene model) with progress display
## - Repeated 10-fold CV (caret)
## - Bootstrap B=500 with progress + ETA
## =========================================================

suppressPackageStartupMessages({
  library(caret)
  library(pROC)
})

## ====== [A] 你只需要改这里：工作目录 ======
WORKDIR <- "D:/360MoveData/Users/57217/Desktop"
setwd(WORKDIR)
cat("Working directory set to:\n", getwd(), "\n\n")

## ====== [B] 文件名（放在 WORKDIR 里） ======
expr_file <- "Combined_Datasets_Matrix_norm.csv"
grp_file  <- "Combined_Datasets_Group.csv"

## ---------- 1) Load expression matrix ----------
cat("[1/6] Loading expression matrix...\n")
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat("  expr dim (genes x samples): ", paste(dim(expr), collapse=" x "), "\n")

# transpose -> rows=samples, cols=genes
X <- as.data.frame(t(expr))
cat("  X dim after transpose (samples x genes): ", paste(dim(X), collapse=" x "), "\n\n")

## ---------- 2) Load group file and align ----------
cat("[2/6] Loading group file and aligning samples...\n")
grp <- read.csv(grp_file, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "group") %in% colnames(grp)))
stopifnot(!anyDuplicated(grp$sample_id))

# keep only samples present in expression matrix
keep <- grp$sample_id %in% rownames(X)
if (!all(keep)) {
  cat("  Dropping ", sum(!keep), " samples in Group not found in Matrix.\n", sep="")
}
grp <- grp[keep, , drop = FALSE]

# align X rows to grp order
X <- X[grp$sample_id, , drop = FALSE]
y <- factor(grp$group, levels = c("AP", "Control"))

cat("  Group counts:\n")
print(table(y))
cat("\n")

## ---------- 3) Subset 4 model genes ----------
cat("[3/6] Subsetting 4 model genes...\n")
genes_model <- c("Vtn", "Anxa3", "Rela", "Iqgap1")
missing_genes <- setdiff(genes_model, colnames(X))
if (length(missing_genes) > 0) {
  stop("Missing genes in expression matrix: ", paste(missing_genes, collapse = ", "))
}
X4 <- X[, genes_model, drop = FALSE]
X4 <- as.data.frame(lapply(X4, as.numeric))
rownames(X4) <- grp$sample_id
cat("  X4 dim (samples x 4 genes): ", paste(dim(X4), collapse=" x "), "\n\n")

## =========================================================
## Part A: Repeated 10-fold CV (caret) + progress
## =========================================================
cat("[4/6] Running repeated 10-fold CV (this may take a bit)...\n")

## （可选）并行加速：如果你机器是多核，强烈建议打开
## Windows 下可用 doParallel
## 如果你不想并行，把下面这段注释掉即可
suppressPackageStartupMessages({
  library(doParallel)
})
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
cat("  Parallel enabled with", n_cores, "cores.\n")

set.seed(123)
ctrl_cv <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter = TRUE          # ✅ 显示每一步训练进度
)

fit_cv <- train(
  x = X4,
  y = y,
  method = "glm",
  family = binomial(),
  metric = "ROC",
  trControl = ctrl_cv
)

stopCluster(cl)
registerDoSEQ()
cat("\n  CV finished.\n")

cv_rocs <- fit_cv$resample$ROC
cat("\n=== Repeated 10-fold CV Summary ===\n")
cat("Mean AUC =", round(mean(cv_rocs), 3),
    " | SD =", round(sd(cv_rocs), 3),
    " | Min =", round(min(cv_rocs), 3),
    " | Max =", round(max(cv_rocs), 3), "\n\n")

## =========================================================
## Part B: Bootstrap B=500 with progress + ETA
## =========================================================
cat("[5/6] Running bootstrap (B=500) with progress...\n")

calc_metrics <- function(y_true, p_hat) {
  y01 <- ifelse(y_true == "AP", 1, 0)
  
  auc <- as.numeric(pROC::auc(pROC::roc(
    y_true, p_hat, levels = c("Control","AP"), direction = "<"
  )))
  
  brier <- mean((p_hat - y01)^2)
  
  eps <- 1e-6
  p_hat2 <- pmin(pmax(p_hat, eps), 1 - eps)
  lp <- qlogis(p_hat2)
  
  fit_int <- suppressWarnings(glm(y01 ~ 1, family = binomial(), offset = lp))
  cal_intercept <- as.numeric(coef(fit_int)[1])
  
  fit_slope <- suppressWarnings(glm(y01 ~ lp, family = binomial()))
  cal_slope <- as.numeric(coef(fit_slope)[2])
  
  list(auc = auc, brier = brier, cal_intercept = cal_intercept, cal_slope = cal_slope)
}

set.seed(123)
B <- 500
n <- nrow(X4)

# full model apparent
fit_full <- glm(y ~ ., data = cbind(y = y, X4), family = binomial())
p_full <- predict(fit_full, type = "response")
m_full <- calc_metrics(y, p_full)

opt_auc   <- numeric(B)
opt_brier <- numeric(B)
opt_int   <- numeric(B)
opt_slope <- numeric(B)

test_auc   <- numeric(B)
test_brier <- numeric(B)
test_int   <- numeric(B)
test_slope <- numeric(B)

t0 <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample.int(n, size = n, replace = TRUE)
  Xb <- X4[idx, , drop = FALSE]
  yb <- y[idx]
  
  fit_b <- suppressWarnings(glm(yb ~ ., data = cbind(yb = yb, Xb), family = binomial()))
  
  # apparent on bootstrap sample
  p_app <- predict(fit_b, type = "response")
  m_app <- calc_metrics(yb, p_app)
  
  # test on original data
  p_test <- predict(fit_b, newdata = X4, type = "response")
  m_t <- calc_metrics(y, p_test)
  
  opt_auc[b]   <- m_app$auc   - m_t$auc
  opt_brier[b] <- m_app$brier - m_t$brier
  opt_int[b]   <- m_app$cal_intercept - m_t$cal_intercept
  opt_slope[b] <- m_app$cal_slope - m_t$cal_slope
  
  test_auc[b]   <- m_t$auc
  test_brier[b] <- m_t$brier
  test_int[b]   <- m_t$cal_intercept
  test_slope[b] <- m_t$cal_slope
  
  ## ✅ 进度显示：每25次报一次
  if (b %% 25 == 0 || b == 1 || b == B) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- elapsed / b
    eta <- rate * (B - b)
    cat(sprintf("  Bootstrap %d/%d (%.1f%%) | elapsed %.1fs | ETA %.1fs\n",
                b, B, 100*b/B, elapsed, eta))
  }
}

# optimism-corrected
auc_corr   <- m_full$auc   - mean(opt_auc)
brier_corr <- m_full$brier - mean(opt_brier)
int_corr   <- m_full$cal_intercept - mean(opt_int)
slope_corr <- m_full$cal_slope - mean(opt_slope)

# percentile 95% CI from test-on-original distribution
ci_auc   <- quantile(test_auc,   c(0.025, 0.975), na.rm = TRUE)
ci_brier <- quantile(test_brier, c(0.025, 0.975), na.rm = TRUE)
ci_int   <- quantile(test_int,   c(0.025, 0.975), na.rm = TRUE)
ci_slope <- quantile(test_slope, c(0.025, 0.975), na.rm = TRUE)

cat("\n=== Bootstrap (B=500) Results ===\n")
cat("[Apparent full-data]\n")
cat("AUC =", round(m_full$auc, 3), "\n")
cat("Brier =", round(m_full$brier, 3), "\n")
cat("Cal intercept =", round(m_full$cal_intercept, 3), "\n")
cat("Cal slope =", round(m_full$cal_slope, 3), "\n\n")

cat("[Test-on-original distribution]\n")
cat("AUC mean =", round(mean(test_auc), 3), " | 95% CI:", round(ci_auc[1],3), "-", round(ci_auc[2],3), "\n")
cat("Brier mean =", round(mean(test_brier), 3), " | 95% CI:", round(ci_brier[1],3), "-", round(ci_brier[2],3), "\n")
cat("Cal intercept mean =", round(mean(test_int), 3), " | 95% CI:", round(ci_int[1],3), "-", round(ci_int[2],3), "\n")
cat("Cal slope mean =", round(mean(test_slope), 3), " | 95% CI:", round(ci_slope[1],3), "-", round(ci_slope[2],3), "\n\n")

cat("[Optimism-corrected]\n")
cat("AUC (corrected) =", round(auc_corr, 3), "\n")
cat("Brier (corrected) =", round(brier_corr, 3), "\n")
cat("Cal intercept (corrected) =", round(int_corr, 3), "\n")
cat("Cal slope (corrected) =", round(slope_corr, 3), "\n")

## ---------- 6) Save results ----------
cat("\n[6/6] Saving results table...\n")
res <- data.frame(
  Metric = c("AUC", "Brier", "Cal_Intercept", "Cal_Slope"),
  Apparent_FullData = c(m_full$auc, m_full$brier, m_full$cal_intercept, m_full$cal_slope),
  Bootstrap_Mean_TestOnOriginal = c(mean(test_auc), mean(test_brier), mean(test_int), mean(test_slope)),
  CI2.5 = c(ci_auc[1], ci_brier[1], ci_int[1], ci_slope[1]),
  CI97.5 = c(ci_auc[2], ci_brier[2], ci_int[2], ci_slope[2]),
  Optimism_Corrected = c(auc_corr, brier_corr, int_corr, slope_corr)
)

write.csv(res, "InternalValidation_4GeneModel_RepeatedCV_B500.csv", row.names = FALSE)
cat("Saved: InternalValidation_4GeneModel_RepeatedCV_B500.csv\n")
cat("\nDONE.\n")



#######################################
# 假设11个基因的名字
genes_all <- c("Pyhin1", "Pah", "Vtn", "Rela", "Mpeg1", "Lcn2", "Car9", "Anxa3", "Actn4", "Iqgap1", "Flna")

# 读取数据并进行基本设置
X <- read.csv("Combined_Datasets_Matrix_norm.csv", row.names = 1)
X <- as.data.frame(t(X))  # 转置：行=样本，列=基因
grp <- read.csv("Combined_Datasets_Group.csv", stringsAsFactors = FALSE)
y <- factor(grp$group, levels = c("AP", "Control"))

# 11个基因的表达数据
X11 <- X[, genes_all]

# 设置bootstrap次数
B <- 500

# 存储每次LASSO选择的基因（选出哪些基因）
selected_genes <- matrix(0, nrow = B, ncol = length(genes_all))
colnames(selected_genes) <- genes_all

# 进行500次bootstrap重采样
# 安装 glmnet 包（如果尚未安装）
install.packages("glmnet")

# 加载 glmnet 包
library(glmnet)


set.seed(123)
for (b in 1:B) {
  cat("Running bootstrap iteration:", b, "/", B, "\n")
  
  # 创建bootstrap样本（随机采样）
  idx <- sample(1:nrow(X11), replace = TRUE)
  Xb <- X11[idx, , drop = FALSE]
  yb <- y[idx]
  
  # 使用 LASSO 进行特征选择
  fit <- cv.glmnet(as.matrix(Xb), yb, alpha = 1, family = "binomial", type.measure = "deviance")
  
  # 选择的非零基因
  selected <- coef(fit, s = "lambda.min")[-1]  # 去掉截距项
  selected_genes[b, ] <- ifelse(selected != 0, 1, 0)  # 标记被选中的基因
}

# 计算每个基因被选中的频率
gene_selection_freq <- colMeans(selected_genes)

# 输出结果
cat("\nGene selection frequency (Bootstrap) for each gene:\n")
print(gene_selection_freq)













