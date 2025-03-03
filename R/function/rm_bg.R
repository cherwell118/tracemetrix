# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'dataName', 'd', 1, 'character','Information of data file.',
  'sampleInfoName', 'a', 1, 'character','Name of Sample file.',
  'dilName', 'l', 1, 'character','Name of dilution file.',
  'method', 'm', 1, 'character','The type of remove data in background: "blank" or "dilution"',
  'fold', 'z', 1, 'numeric','Fold number',
  'r2', 'c', 1, 'numeric','The number of r2 cutoff'
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
check_required_param("method", opt$method, "Method (blank or dilution) is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
method <- opt$method

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

sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
sampleInfo <- read_csv(sampleInfoFile, col_types = cols()) %>% as.data.frame()
sampleInfo <- sampleInfo[!apply(sampleInfo, 1, function(x) all(is.na(x))), ]
sampleInfo <- sampleInfo[, !apply(sampleInfo, 2, function(x) all(is.na(x)))]
if (nrow(sampleInfo) == 0) {
  stop("The sample information file is empty!")
}

# 复制样本信息文件到输出目录
if (file.exists(paste0(outputDir, "/", sampleInfoName))) {
  file.remove(paste0(outputDir, "/", sampleInfoName))
}
file.copy(from = sampleInfoFile, to = paste0(outputDir, "/", sampleInfoName))

# 背景扣除
if (method == "blank") {
  fold <- ifelse(is.null(opt$fold), 3, as.numeric(opt$fold))
  if (!("Blank" %in% unique(sampleInfo$classes))) {
    stop("'Blank' is not contained in column 'classes' of sample info file!")
  }
  classes <- ifelse("QC" %in% unique(sampleInfo$classes), c("QC", "Sample"), c("Sample"))
  cal_index <- lapply(classes, function(x) get_index(colnames(data), sampleInfo, group = x, by = "classes")) %>% unlist
  blank_index <- get_index(colnames(data), sampleInfo, group = "Blank", by = "classes")
  
  data_filted <- bk_deduct_oneBlank(data, blank_index, cal_index, fold = fold)
} else if (method == "dilution") {
  r2 <- ifelse(is.null(opt$r2), 0.8, as.numeric(opt$r2))
  dilName <- opt$dilName
  dilFile <- paste0(inputDir, "/", dilName)
  if (!file.exists(dilFile)) {
    stop(paste0(dilFile, " does not exist!"))
  }
  dil_data <- read_csv(dilFile, col_types = cols()) %>% arrange(order) %>% as.data.frame()
  dil_data <- dil_data[!apply(dil_data, 1, function(x) all(is.na(x))), ]
  dil_data <- dil_data[, !apply(dil_data, 2, function(x) all(is.na(x)))]
  
  diluted_data <- data[, dil_data$sample.names]
  diluted_data[apply(diluted_data, 1, function(x) all(is.na(x))), ] <- 0
  line_data <- apply(diluted_data, 1, function(x) linear(x, dilution = dil_data$dilution))
  line_data <- t(line_data)
  colnames(line_data) <- c("k", "b", "r2")
  
  data <- cbind(data, line_data)
  data_filted <- data[data$k > 0 & data$r2 >= r2, ]
  cat(paste0("The feature size after filtering by k and r2 is ", nrow(data_filted), ".\n"))
  cat(paste0(nrow(data) - nrow(data_filted), " features have been removed.\n"))
} else {
  stop("The parameter 'method' must be 'blank' or 'dilution'.")
}

# RSD 计算
if (is.null(dim(data_filted))) {
  stop("Data after filtering is empty.")
}

sample_index <- get_index(colnames(data_filted), sampleInfo, group = "Sample", by = "classes")
sample <- data_filted[, sample_index]
sample.rsd <- apply(sample, 1, rsdFun)
data_filted$sample.rsd <- sample.rsd

if ("QC" %in% unique(sampleInfo$classes)) {
  qc_index <- get_index(colnames(data_filted), sampleInfo, group = "QC", by = "classes")
  qc <- data_filted[, qc_index]
  QC.rsd <- apply(qc, 1, rsdFun)
  data_filted$QC.rsd <- QC.rsd
}

# 输出结果
write.csv(data_filted, paste0(outputDir, "/data.csv"), row.names = FALSE)
plotDir <- paste0(outputDir, "/plotData")
if (!dir.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
suppressWarnings(qcplot(inputDir = outputDir, outputDir = plotDir, dataName = "data.csv", sampleInfoName = sampleInfoName))