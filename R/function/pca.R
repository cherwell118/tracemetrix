# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, "logical",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'sampleInfoName', 'a', 1, "character",'Name of Sample file.',
  'scaleC', 'S', 1, 'character','One of none, center, standard and pareto',
  'group', 'g', 1, 'character','Name of column which used as grouping information'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

# 帮助信息
if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
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
check_required_param("scaleC", opt$scaleC, "Scaling method is required.")
check_required_param("group", opt$group, "Grouping column name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
scaleC <- opt$scaleC
group <- opt$group

# 创建输出目录
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# 加载数据
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]
if (nrow(data) == 0) {
  stop("The data file is empty!")
}
rownames(data) <- data$compound.id

sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
sampleInfo <- read_csv(sampleInfoFile, col_types = cols()) %>% as.data.frame()
sampleInfo <- sampleInfo[!apply(sampleInfo, 1, function(x) all(is.na(x))), ]
sampleInfo <- sampleInfo[, !apply(sampleInfo, 2, function(x) all(is.na(x)))]
if (nrow(sampleInfo) == 0) {
  stop("The sample information file is empty!")
}
sampleInfo <- sampleInfo[sampleInfo$sample.names %in% colnames(data), ]

# PCA 函数
perform_pca <- function(data, sampleInfo, group, scaleC, outputDir) {
  groups <- unique(sampleInfo[[group]])
  cal_index <- lapply(groups, function(x) get_index(colnames(data), sampleInfo, group = x, by = group)) %>% unlist
  cal_index <- sort(cal_index)
  index_data <- data.frame(groups = names(cal_index), index = cal_index)

  pca_dat <- data[, index_data$index]

  # 数据标准化
  if (scaleC == "none") {
    pca_dat <- t(pca_dat)
  } else if (scaleC == "center") {
    pca_dat <- scale(t(pca_dat), center = TRUE, scale = FALSE)
  } else if (scaleC == "standard") {
    pca_dat <- scale(t(pca_dat), center = TRUE, scale = TRUE)
  } else if (scaleC == "pareto") {
    pca_dat <- scale(t(pca_dat), center = TRUE, scale = FALSE)
    pca_dat <- apply(pca_dat, 2, function(colData) { colData / sqrt(sd(colData)) })
  }

  write.csv(pca_dat, paste0(outputDir, "/data_", scaleC, ".csv"))

  pca_object <- prcomp(pca_dat, center = FALSE, scale = FALSE)
  write.csv(pca_object$x, paste0(outputDir, "/PCA_components_", scaleC, ".csv"))
  contr <- summary(pca_object)$importance
  write.csv(contr, paste0(outputDir, "/PCA_contr.csv"))

  pca_loading <- pca_object$rotation
  write.csv(pca_loading, paste0(outputDir, "/PCA_loading_", scaleC, ".csv"))

  return(list(pca_object = pca_object, contr = contr, index_data = index_data))
}

# 执行 PCA
pca_results <- perform_pca(data, sampleInfo, group, scaleC, outputDir)
pca_object <- pca_results$pca_object
contr <- pca_results$contr
index_data <- pca_results$index_data

# 可视化 PCA 结果
plot_pca <- function(pca_object, contr, index_data, scaleC, outputDir) {
  pcn <- 8
  pv <- contr[2, ]
  if (length(pv) > pcn) {
    pv <- pv[1:pcn]
  }
  plot_data <- pca_object$x[, 1:pcn]
  var_labels <- paste0(colnames(plot_data)[1:pcn], " (", round(pv * 100, 2), "%)")

  my_colors <- hcl.colors(nlevels(as.factor(index_data$groups)), palette = "viridis")

  pca_df <- data.frame(pca_object$x[, 1:2], groups = as.character(index_data$groups))

  pc1_lab <- paste0("PC1 (", round(contr[2, 1] * 100, 2), "%)")
  pc2_lab <- paste0("PC2 (", round(contr[2, 2] * 100, 2), "%)")

  pdf_file <- paste0(outputDir, "/PCA_", scaleC, ".pdf")
  pdf(pdf_file)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = groups)) +
    scale_colour_manual(values = my_colors) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15)) +
    geom_point() +
    xlab(pc1_lab) +
    ylab(pc2_lab)
  print(p)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = groups)) +
    scale_colour_viridis_d() +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15)) +
    geom_point() +
    xlab(pc1_lab) +
    ylab(pc2_lab) +
    stat_ellipse(show.legend = FALSE, level = 0.95)
  print(p)
  dev.off()

  write.csv(t(contr), paste0(outputDir, "/principle.Variance_", scaleC, ".csv"))

  # 碎石图
  pv <- contr[2, ]
  pv <- data.frame(pv = pv, x = 1:length(pv))
  if (nrow(pv) > 20) {
    pv <- pv[1:20, ]
  }

  p <- ggplot(pv, aes(x = x, y = pv)) +
    theme_bw() +
    geom_point() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15)) +
    geom_line(col = "grey", lwd = 1) +
    geom_point(pch = 16, lwd = 2, col = 1) +
    ylab("Ratio of Variances") +
    xlab("PCn")

  pdf_file <- paste0(outputDir, "/screeplot.Ratio_of_Variances_", scaleC, ".pdf")
  pdf(pdf_file)
  print(p)
  dev.off()

  pdf_file <- paste0(outputDir, "/screeplot.Variances_", scaleC, ".pdf")
  pdf(pdf_file)
  screeplot(pca_object, type = "lines", main = "Scree Plot")
  dev.off()
}

# 绘制 PCA 图
plot_pca(pca_object, contr, index_data, scaleC, outputDir)