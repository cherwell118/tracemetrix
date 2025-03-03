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
  'comparisionName', 'c', 1, 'character','The comparison data file.',
  'method', 'm', 1, 'character','Method: t.test, wilcox.test, or auto.',
  'group', 'g', 1, 'character','Name of column used for grouping.',
  'var.equal', 'v', 1, 'logical','Logical variable for t.test indicating equal variances.',
  'alternative', 't', 1, 'character','Alternative hypothesis: "two.sided", "greater", or "less".'
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
check_required_param("comparisionName", opt$comparisionName, "Comparison file name is required.")
check_required_param("method", opt$method, "Method (t.test, wilcox.test, or auto) is required.")
check_required_param("group", opt$group, "Group column name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
comparisionName <- opt$comparisionName
method <- opt$method
group <- opt$group
var.equal <- ifelse(is.null(opt$var.equal), TRUE, as.logical(opt$var.equal))
alternative <- ifelse(is.null(opt$alternative), "two.sided", opt$alternative)

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
sampleInfo <- sampleInfo[sampleInfo$sample.names %in% colnames(data), ]

comparisionFile <- paste0(inputDir, "/", comparisionName)
comparision <- read_csv(comparisionFile, col_types = cols()) %>% as.data.frame()
comparision <- comparision[!apply(comparision, 1, function(x) all(is.na(x))), ]
comparision <- comparision[, !apply(comparision, 2, function(x) all(is.na(x)))]
if (!("paired" %in% colnames(comparision))) {
  comparision$paired <- rep(FALSE, nrow(comparision))
}

# 比较函数
compare <- function(x, comparision, sampleInfo, data, method, var.equal, alternative) {
  groupA <- comparision[x, 1]
  groupB <- comparision[x, 2]
  paired <- comparision$paired[x] %>% as.logical
  
  infoA <- sampleInfo[sampleInfo[[group]] == groupA, ]
  infoB <- sampleInfo[sampleInfo[[group]] == groupB, ]
  
  if (paired) {
    if (!identical(sort(infoA$label), sort(infoB$label))) {
      stop("There is a problem with the paired information.")
    }
    nameA <- infoA$sample.names[order(infoA$label)]
    nameB <- infoB$sample.names[order(infoB$label)]
  } else {
    nameA <- infoA$sample.names
    nameB <- infoB$sample.names
  }
  
  dataA <- data[, nameA] %>% as.matrix() %>% as.numeric() %>% matrix(ncol = length(nameA))
  dataB <- data[, nameB] %>% as.matrix() %>% as.numeric() %>% matrix(ncol = length(nameB))
  
  rownames(dataA) <- rownames(dataB) <- data$compound.id
  colnames(dataA) <- nameA
  colnames(dataB) <- nameB
  
  if (method == "auto") {
    pvalue <- get_p(dataA, dataB, method = "auto", paired = paired, alternative = alternative) %>% as.numeric()
  } else if (method == "t.test") {
    pvalue <- get_p(dataA, dataB, method = "t.test", var.equal = var.equal, paired = paired, alternative = alternative) %>% as.numeric()
  } else if (method == "wilcox.test") {
    pvalue <- get_p(dataA, dataB, method = "wilcox.test", paired = paired, alternative = alternative) %>% as.numeric()
  }
  
  MeanA <- apply(dataA, 1, mean)
  MeanB <- apply(dataB, 1, mean)
  meanDF <- data.frame(MeanA = MeanA, MeanB = MeanB)
  colnames(meanDF) <- c(paste0("mean_", groupA), paste0("mean_", groupB))
  
  fdr <- p.adjust(pvalue, "BH") %>% as.numeric()
  bon <- p.adjust(pvalue, "bonferroni") %>% as.numeric()
  log2FC <- log2(MeanA / MeanB)
  
  results <- data.frame(log2FC, pvalue, fdr, bon)
  colnames(results) <- c("log2FC", "pvalue", "fdr", "bonferroni")
  results <- cbind(dataA, dataB, meanDF, results)
  
  return(results)
}

# 执行比较
sink(paste0(outputDir, "/cal_pvalue_formula.txt"))
compare_results <- lapply(1:nrow(comparision), function(x) {
  cat(paste0(comparision[x, 1], "_vs_", comparision[x, 2], " ", Sys.time(), "\n"))
  compare(x, comparision, sampleInfo, data, method, var.equal, alternative)
})
sink()

# 输出结果
for (i in 1:nrow(comparision)) {
  results <- compare_results[[i]]
  results <- cbind(Compound = rownames(results), compare_results[[i]])
  write.csv(results, paste0(outputDir, "/results_", comparision[i, 1], "_vs_", comparision[i, 2], ".csv"), row.names = FALSE)
}