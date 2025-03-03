# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(mice)))
suppressPackageStartupMessages(suppressWarnings(library(missForest)))
suppressPackageStartupMessages(suppressWarnings(library(imputation)))
suppressPackageStartupMessages(suppressWarnings(library(pcaMethods)))
suppressPackageStartupMessages(suppressWarnings(library(impute)))
suppressPackageStartupMessages(suppressWarnings(library(SeqKnn)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'dataName', 'd', 1, 'character','Information of data file.',
  'sampleInfoName', 'a', 1, 'character','Name of Sample file.',
  'naRatio', 'r', 1, 'numeric','NA Ratio cutoff for feature removal',
  'kNum', 'k', 1, 'character','Parameter in svd, seqknn, and knn method.',
  'scaleMethod', 's', 1, 'character','Scale method in bpca or ppca.',
  'method', 'm', 1, 'character','Method for missing value process.'
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
check_required_param("method", opt$method, "Method for missing value process is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
naRatio <- ifelse(is.null(opt$naRatio), 0.2, as.numeric(opt$naRatio))
method <- opt$method
kNum <- ifelse(is.null(opt$kNum), 5, as.numeric(opt$kNum))
scaleMethod <- ifelse(is.null(opt$scaleMethod), "none", opt$scaleMethod)

# 创建输出目录
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# 复制样本信息文件到输出目录
sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
outputSampleInfoFile <- paste0(outputDir, "/", sampleInfoName)

if (file.exists(outputSampleInfoFile)) {
  file.remove(outputSampleInfoFile)
}
file.copy(from = sampleInfoFile, to = outputSampleInfoFile)

# 加载数据
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]

sampleInfo <- read_csv(sampleInfoFile, col_types = cols()) %>% as.data.frame()
sampleInfo <- sampleInfo[!apply(sampleInfo, 1, function(x) all(is.na(x))), ]
sampleInfo <- sampleInfo[, !apply(sampleInfo, 2, function(x) all(is.na(x)))]

# 特征筛选：根据NA比例删除特征
feature_rm <- function(data, ratio, rm.type = "NA") {
  na_ratio <- apply(data, 1, function(x) sum(is.na(x)) / length(x))
  keep_index <- which(na_ratio <= ratio)
  return(list(row_index = keep_index, na_ratio = na_ratio))
}

data_rm <- feature_rm(data, ratio = naRatio, rm.type = "NA")
data_filted <- data[data_rm$row_index, ]

if (is.null(dim(data_filted))) {
  stop("Data after filtering is empty.")
}

# 缺失值处理
process_missing_values <- function(data, method, kNum, scaleMethod) {
  cat(paste0(Sys.time(), " Start processing missing values using ", method, "\n"))
  switch(method,
         randomForest = {
           set.seed(123)
           data <- t(data)
           data <- missForest(data, ntree = 100)[[1]] %>% t
         },
         SVD = {
           data <- SVDImpute(data, k = kNum, num.iters = 10, verbose = FALSE)
           data <- data$x
         },
         PPCA = {
           data_imp <- pca(data, method = "ppca", nPcs = 4, scale = scaleMethod, center = TRUE)
           data <- completeObs(data_imp) %>% as.data.frame()
         },
         BPCA = {
           data_imp <- pca(data, method = "bpca", nPcs = 4, scale = scaleMethod, center = TRUE)
           data <- completeObs(data_imp) %>% as.data.frame()
         },
         KNN = {
           data_imp <- impute.knn(t(data), k = kNum, colmax = 0.8, rowmax = 0.8, rng.seed = 123)
           data <- t(data_imp$data) %>% as.data.frame()
         },
         SKNN = {
           data <- SeqKNN(data, k = kNum)
         },
         mean = {
           data <- apply(data, 1, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)) %>% t
         },
         median = {
           data <- apply(data, 1, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)) %>% t
         },
         minimum = {
           data <- apply(data, 1, function(x) ifelse(is.na(x), min(x, na.rm = TRUE), x)) %>% t
         },
         halfMinimum = {
           data <- apply(data, 1, function(x) ifelse(is.na(x), min(x, na.rm = TRUE) / 2, x)) %>% t
         },
         zero = {
           data <- apply(data, 1, function(x) ifelse(is.na(x), 0, x)) %>% t
         })
  cat(paste0(Sys.time(), " End of processing missing values using ", method, "\n"))
  return(data)
}

data_mv <- process_missing_values(data_filted, method, kNum, scaleMethod)

# 保存处理后的数据
data_info <- data_filted[, which(colnames(data_filted) %in% c("compound.id", "mz", "rt")), drop = FALSE]
data_final <- cbind(data_info, data_mv)
write.csv(data_final, paste0(outputDir, "/data.csv"), row.names = FALSE)

# 保存样本信息
sampleInfo <- sampleInfo[sampleInfo$classes %in% c("Sample", "QC"), ]
write.csv(sampleInfo, paste0(outputDir, "/sample_info.csv"), row.names = FALSE)

# RSD计算
if ("Sample" %in% unique(sampleInfo$classes)) {
  sample_index <- get_index(colnames(data_final), sampleInfo, group = "Sample", by = "classes")
  sample <- data_final[, sample_index]
  sample.rsd <- apply(sample, 1, rsdFun)
  data_final$sample.rsd <- sample.rsd
}

if ("QC" %in% unique(sampleInfo$classes)) {
  qc_index <- get_index(colnames(data_final), sampleInfo, group = "QC", by = "classes")
  qc <- data_final[, qc_index]
  QC.rsd <- apply(qc, 1, rsdFun)
  data_final$QC.rsd <- QC.rsd
}

write.csv(data_final, paste0(outputDir, "/data_with_rsd.csv"), row.names = FALSE)

# 可视化
plotDir <- paste0(outputDir, "/plotData")
if (!dir.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

suppressWarnings(qcplot(inputDir = outputDir, outputDir = plotDir, dataName = "data.csv", sampleInfoName = "sample_info.csv"))