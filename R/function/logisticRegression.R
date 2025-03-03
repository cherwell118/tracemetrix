# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(car)))
suppressPackageStartupMessages(suppressWarnings(library(visreg)))
suppressPackageStartupMessages(suppressWarnings(library(survival)))
suppressPackageStartupMessages(suppressWarnings(library(MASS)))
suppressPackageStartupMessages(suppressWarnings(library(nnet)))
suppressPackageStartupMessages(suppressWarnings(library(xlsx)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'method', 'm', 1, "character",'One of binary, polytomous and ordinal',
  'y', 'y', 1, "character",'Name of dependent variable.',
  'x', 'x', 1, "character",'Name of independent variable(s).',
  'correction', 'R', 1, "character",'Name(s) of variable(s) for correction, comma separated.',
  'categorical', 'C', 1, "character",'Names of categorical variable(s), comma separated.'
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
check_required_param("y", opt$y, "Dependent variable name is required.")
check_required_param("x", opt$x, "Independent variable(s) name is required.")
check_required_param("method", opt$method, "Method (binary, polytomous, ordinal) is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
y <- opt$y
x <- opt$x
method <- opt$method
correction <- opt$correction
categorical <- opt$categorical

# 创建输出目录
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# 加载数据
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]

# 数据预处理
prepare_data <- function(data, x, categorical) {
  independent <- unlist(strsplit(x, split = ","))
  if (!is.null(categorical)) {
    categorical <- unlist(strsplit(categorical, split = ","))
    for (variable in categorical) {
      if (!is.factor(data[[variable]])) {
        data[[variable]] <- as.factor(data[[variable]])
      }
    }
  }
  return(list(data = data, independent = independent))
}

data_prepared <- prepare_data(data, x, categorical)
data <- data_prepared$data
independent <- data_prepared$independent

# 生成公式
generate_formula <- function(y, independent, correction = NULL) {
  if (is.null(correction)) {
    formula <- as.formula(paste(y, "~", paste(independent, collapse = "+")))
  } else {
    correction <- unlist(strsplit(correction, split = ","))
    formula <- as.formula(paste(y, "~", paste(c(correction, independent), collapse = "+")))
  }
  return(formula)
}

fmla <- generate_formula(y, independent, correction)

# 逻辑回归分析
logistic_regression <- function(data, formula, method, outputDir) {
  if (method == "binary") {
    fit <- glm(formula, family = binomial(link = "logit"), data = data)
  } else if (method == "polytomous") {
    fit <- multinom(formula, data = data)
  } else if (method == "ordinal") {
    fit <- polr(formula, data = data, Hess = TRUE, method = "logistic")
  } else {
    stop("Invalid method. Choose from 'binary', 'polytomous', or 'ordinal'.")
  }
  
  # 保存摘要信息
  sink_file <- paste0(outputDir, "/summary.txt")
  sink(sink_file)
  print(summary(fit))
  sink()
  
  # 提取结果
  if (method == "binary") {
    coefficients <- coef(summary(fit))
    write.csv(coefficients, paste0(outputDir, "/coefficient_table.csv"), row.names = TRUE)
  } else if (method == "polytomous") {
    coefficients <- coef(fit)
    write.csv(coefficients, paste0(outputDir, "/coefficient_table.csv"), row.names = TRUE)
  } else if (method == "ordinal") {
    coefficients <- coef(summary(fit))
    write.csv(coefficients, paste0(outputDir, "/coefficient_table.csv"), row.names = TRUE)
  }
  
  return(fit)
}

fit <- logistic_regression(data, fmla, method, outputDir)

# 可视化
plot_results <- function(fit, data, method, outputDir) {
  if (method == "binary") {
    pdf_file <- paste0(outputDir, "/regressionPlot.pdf")
    pdf(pdf_file)
    ggplot(data, aes(x = data[, independent[1]], y = data[, y])) +
      geom_point() +
      stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
      ylab(y) +
      xlab(independent[1]) +
      theme_bw()
    dev.off()
  } else if (method == "polytomous") {
    pdf_file <- paste0(outputDir, "/regressionPlot.pdf")
    pdf(pdf_file)
    visreg::visreg(fit)
    dev.off()
  } else if (method == "ordinal") {
    pdf_file <- paste0(outputDir, "/regressionPlot.pdf")
    pdf(pdf_file)
    visreg::visreg(fit)
    dev.off()
  }
}

plot_results(fit, data, method, outputDir)

# 异常点检测（仅限二分类逻辑回归）
if (method == "binary") {
  pdf_file <- paste0(outputDir, "/influencePlot.pdf")
  pdf(pdf_file)
  car::infIndexPlot(fit)
  points <- car::influencePlot(fit)
  dev.off()
  
  write.csv(points, paste0(outputDir, "/influencePoint.csv"))
  write.csv(car::infIndex(fit), paste0(outputDir, "/influencePoint.all.csv"))
}

# 过度离势检测（仅限二分类逻辑回归）
if (method == "binary") {
  p <- pchisq(summary(fit)$dispersion * fit$df.residual, fit$df.residual, lower = FALSE)
  sink_file <- paste0(outputDir, "/overdispersion_report.txt")
  sink(sink_file)
  if (p >= 0.05) {
    cat(paste0("p = ", signif(p, digits = 4), " >= 0.05, there is no overdispersion\n"))
  } else {
    cat(paste0("p = ", signif(p, digits = 4), " < 0.05, there is overdispersion\n"))
  }
  sink()
}