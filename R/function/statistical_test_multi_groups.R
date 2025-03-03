# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, "logical",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'sampleInfoName', 'a', 1, "character",'Name of Sample file.',
  'method', 'm', 1, 'character','Method: standar.anova, Welch.anova, or auto',
  'group', 'g', 1, 'character','Name of column used for grouping information'
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
check_required_param("method", opt$method, "Method (standar.anova, Welch.anova, or auto) is required.")
check_required_param("group", opt$group, "Group column name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
method <- opt$method
group <- opt$group

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

# 获取分组信息
groups <- unique(sampleInfo[[group]])
cal_index <- lapply(groups, function(x) get_index(colnames(data), sampleInfo, group = x, by = group)) %>% unlist
index <- data.frame(groups = names(cal_index), index = cal_index)
index <- index[order(index$index), ]
cal_data <- data[, index$index]
cal_data <- cbind(index = c(1:nrow(cal_data)), cal_data)

# 计算均值和标准差
multiMean <- function(data, groups) {
  means <- tapply(data, INDEX = groups, FUN = mean)
  return(means)
}

multiSd <- function(data, groups) {
  sds <- tapply(data, INDEX = groups, FUN = sd)
  return(sds)
}

means <- apply(cal_data, 1, function(x) {
  idx <- as.numeric(x)[1]
  y <- as.numeric(x)[-1]
  results <- multiMean(data = y, groups = index$groups)
  return(results)
})

sd <- apply(cal_data, 1, function(x) {
  idx <- as.numeric(x)[1]
  y <- as.numeric(x)[-1]
  results <- multiSd(data = y, groups = index$groups)
  return(results)
})

means <- t(means)
colnames(means) <- paste0(colnames(means), "_mean")

# 执行多组样本的统计检验
multiTest <- function(data, groups, method) {
  if (method == "standar.anova") {
    aov_result <- aov(data ~ groups)
    summary_result <- summary(aov_result)
    p_value <- summary_result[[1]]$`Pr(>F)`[1]
    method_name <- "standar anova"
  } else if (method == "Welch.anova") {
    aov_result <- oneway.test(data ~ groups, var.equal = FALSE)
    summary_result <- summary(aov_result)
    p_value <- summary_result[[1]]$`Pr(>F)`[1]
    method_name <- "Welch anova"
  } else {
    aov_result <- aov(data ~ groups)
    summary_result <- summary(aov_result)
    p_value <- summary_result[[1]]$`Pr(>F)`[1]
    method_name <- "auto"
  }
  return(list(p = p_value, method = method_name, table = summary_result[[1]]))
}

sink(paste0(outputDir, "/cal_pvalue_methods.txt"))
p <- apply(cal_data, 1, function(x) {
  idx <- as.numeric(x)[1]
  cat(paste0(rownames(cal_data)[idx], ":\n"))
  y <- as.numeric(x)[-1]
  results <- multiTest(data = y, groups = index$groups, method = method)
  anova_p <- results$p
  if (results$method == "standar anova") {
    p <- results$table[,"p adj"]
    method <- 1
  } else {
    p <- results$table[,"p"]
    method <- 2
  }
  if (length(p) == 1) {
    names(p) <- rownames(results$table)
  }
  p <- c(method = method, anova_p = anova_p, p)
  return(p)
})
sink()

p %<>% t(.) %>% as.data.frame
p$method <- ifelse(p$method == 1, "standar anova", "Welch anova")

info <- p[, 1:2]
p <- p[, -c(1:2), drop = FALSE]

info$anova_p_fdr <- p.adjust(info$anova_p, "BH")
info$anova_p_bon <- p.adjust(info$anova_p, "bonferroni")

col_index <- do.call(paste, c(expand.grid(colnames(p), c("pvalue", "fdr", "bon")) %>% arrange(Var1), sep = "_"))

fdr <- apply(p, 2, function(x) { p.adjust(x, "BH") })
colnames(fdr) <- paste0(colnames(fdr), "_fdr")

bon <- apply(p, 2, function(x) { p.adjust(x, "bonferroni") })
colnames(bon) <- paste0(colnames(bon), "_bon")

colnames(p) <- paste0(colnames(p), "_pvalue")

compare_results <- cbind(p, fdr, bon)
compare_results <- compare_results[, col_index]

# 输出结果
data <- cbind(data, means, info, compare_results)
write.csv(data, paste0(outputDir, "/multiple_test_results.csv"), row.names = FALSE)