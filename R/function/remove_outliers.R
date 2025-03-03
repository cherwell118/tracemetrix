# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(car)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'sampleInfoName', 'a', 1, "character",'Name of Sample file.',
  'method', 'm', 1, 'character','One of none, centering, autoscaling and paretoscaling',
  'CI', 'c', 1, 'character', 'Confidence interval, default is 95%'
), byrow = TRUE, ncol = 5)

opt = getopt(commands)

# 帮助信息
if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status = 1)
}

# 参数检查
check_required_param <- function(param_name, param_value, error_message) {
  if (is.null(param_value)) {
    stop(paste0("Please set the parameter '", param_name, "': ", error_message))
  }
}

check_required_param("inputDir", opt$inputDir, "Input directory is required.")
check_required_param("outputDir", opt$outputDir, "Output directory is required.")
check_required_param("dataName", opt$dataName, "Data file name is required.")
check_required_param("sampleInfoName", opt$sampleInfoName, "Sample information file name is required.")
check_required_param("method", opt$method, "Scaling method is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
method <- opt$method
CI <- ifelse(is.null(opt$CI), 0.95, as.numeric(opt$CI))

# 创建输出目录
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# 加载数据
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
if (nrow(data) == 0) {
  stop("The data file is empty!")
}
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]
rownames(data) <- data$Compound

sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
sampleInfo <- read_csv(sampleInfoFile, col_types = cols()) %>% as.data.frame()
if (nrow(sampleInfo) == 0) {
  stop("The sample information file is empty!")
}
sampleInfo <- sampleInfo[!apply(sampleInfo, 1, function(x) all(is.na(x))), ]
sampleInfo <- sampleInfo[, !apply(sampleInfo, 2, function(x) all(is.na(x)))]

# PCA 数据准备
classes <- ifelse("QC" %in% unique(sampleInfo$classes), c("QC", "Sample"), c("Sample"))
cal_index <- lapply(classes, function(x) get_index(colnames(data), sampleInfo, group = x, by = "classes")) %>% unlist
cal_index <- sort(cal_index)
index_data <- data.frame(groups = names(cal_index), index = cal_index)
pca_dat <- data[, index_data$index]

# 数据标准化
scale_data <- function(data, method) {
  if (method == "none") {
    return(data)
  } else if (method == "centering") {
    return(scale(t(data), center = TRUE, scale = FALSE))
  } else if (method == "autoscaling") {
    return(scale(t(data), center = TRUE, scale = TRUE))
  } else if (method == "paretoscaling") {
    scaled_data <- scale(t(data), center = TRUE, scale = FALSE)
    return(apply(scaled_data, 2, function(colData) { colData / sqrt(sd(colData)) }))
  }
}

pca_dat <- scale_data(pca_dat, method)
write.csv(t(pca_dat), paste0(outputDir, "/data_", method, ".csv"))

# PCA 分析
pca_object <- prcomp(pca_dat, center = FALSE, scale = FALSE)
scores <- pca_object$x[, 1:8]
write.csv(scores, paste0(outputDir, "/pca_scores.csv"))

# Hotelling's T2 椭圆计算
HotE <- function(x, y, len = 200, alfa = CI) {
  N <- length(x)
  alfa <- as.numeric(alfa)
  A <- 2
  mypi <- seq(0, 2 * pi, length = len)
  r1 <- sqrt(var(x) * qf(alfa, 2, N - 2) * (2 * (N^2 - 1) / (N * (N - 2))))
  r2 <- sqrt(var(y) * qf(alfa, 2, N - 2) * (2 * (N^2 - 1) / (N * (N - 2))))
  cbind(r1 * cos(mypi) + mean(x), r2 * sin(mypi) + mean(y))
}

# 异常值检测
detect_outliers <- function(scores, pc1, pc2, CI) {
  HotEllipse <- abs(cbind(HotE(scores[, pc1], scores[, pc2])))
  outliers <- as.numeric()
  for (i in 1:nrow(scores)) {
    sample <- abs(scores[i, c(pc1, pc2)])
    out_pc1 <- which(HotEllipse[, 1] < sample[1])
    outlier <- any(HotEllipse[out_pc1, 2] < sample[2]) * 1
    outliers <- c(outliers, outlier)
  }
  return(outliers)
}

# 计算所有 PC 组合的异常值
pc_pairs <- list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))
outliers_results <- lapply(pc_pairs, function(pair) {
  outliers <- detect_outliers(scores, pair[1], pair[2], CI)
  return(data.frame(PC1 = scores[, pair[1]], PC2 = scores[, pair[2]], outliers = outliers == 1))
})

# 保存异常值结果
outlier_names <- lapply(outliers_results, function(result) {
  rownames(result)[result$outliers == TRUE]
})

names(outlier_names) <- c("PC1_PC2", "PC1_PC3", "PC1_PC4", "PC2_PC3", "PC2_PC4", "PC3_PC4")
capture.output(outlier_names, file = paste0(outputDir, "/outliers.txt"))

# 绘制 PCA 得分图
plot_pca_scores <- function(scores, pc_pairs, CI, outputDir) {
  for (i in seq_along(pc_pairs)) {
    pair <- pc_pairs[[i]]
    pc1_lab <- paste0("PC", pair[1], " (", round(summary(pca_object)$importance[2, pair[1]] * 100, 2), "%)")
    pc2_lab <- paste0("PC", pair[2], " (", round(summary(pca_object)$importance[2, pair[2]] * 100, 2), "%)")
    
    pdf_file <- paste0(outputDir, "/outliers_by_PC", pair[1], "_and_PC", pair[2], ".pdf")
    pdf(pdf_file, height = 6, width = 8.5)
    
    p <- ggplot() +
      geom_polygon(data = outliers_results[[i]], aes(PC1, PC2), color = "lightcyan", fill = "lightcyan") +
      theme_bw() +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            strip.background = element_rect(fill = "#0099B47F"),
            strip.text = element_text(color = "white", size = 15)) +
      xlab(pc1_lab) + ylab(pc2_lab) +
      geom_hline(yintercept = 0, lty = 2, col = "red", lwd = 1) +
      geom_vline(xintercept = 0, lty = 2, col = "red", lwd = 1) +
      geom_point(data = outliers_results[[i]], aes(x = PC1, y = PC2, color = outliers), size = 3) +
      scale_color_manual(values = c("darkturquoise", "#3D3D3D"))
    
    print(p)
    dev.off()
  }
}

plot_pca_scores(scores, pc_pairs, CI, outputDir)

# 保存数据和样本信息
write.csv(data, paste0(outputDir, "/data.csv"))
write.csv(sampleInfo, paste0(outputDir, "/sample_info.csv"))