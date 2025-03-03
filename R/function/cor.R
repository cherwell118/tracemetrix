suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "logical",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName1','d',1,"character",'Name of data file 1',
  'dataName2','D',1,"character",'Name of data file 2',
  'method','m','1','character','Method used in Correlation analysis.One of "pearson", "kendall", "spearman"',
  'alternative','t',1,'character','One of "two.sided(default)" , "greater" or "less".'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
}

# Extract options
inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName1 <- opt$dataName1
dataName2 <- opt$dataName2
method <- opt$method
alternative <- opt$alternative

# Default values
if (is.null(method)) {
  stop("Please choose a method")
}
if (is.null(alternative)) {
  alternative <- "two.sided"
}

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Function to load and clean data
load_and_clean_data <- function(dataFile) {
  data <- read_csv(dataFile,
                   col_types =(),
 cols                   progress = FALSE) %>% as.data.frame()
  data <- data[!apply(data, 1, function(x) { all(is.na(x)) }), ]
  data <- data[, !apply(data, 2, function(x) { all(is.na(x)) })]
  if (nrow(data) == 0) {
    stop("The data file is empty!")
  }
  data <- data[order(data$name), ]
  rownames(data) <- data$name
  data <- data[, -1]
  return(data)
}

# Load data files
dataFile <- paste0(inputDir, "/", dataName1)
data <- load_and_clean_data(dataFile)

if (!is.null(dataName2)) {
  dataFile2 <- paste0(inputDir, "/", dataName2)
  data2 <- load_and_clean_data(dataFile2)
  data <- data[(rownames(data) %in% rownames(data2)), ]
  data2 <- data2[(rownames(data2) %in% rownames(data)), ]
}

# Function for correlation analysis
correlation_analysis <- function(data1, data2 = NULL, method, alternative) {
  cat(paste0(Sys.time(), " Start of cor.test.\n"))
  if (is.null(data2)) {
    results <- cor.test_matrix(data1, alternative = alternative, method = method)
  } else {
    results <- cor.test_matrix(data1, data2, alternative = alternative, method = method)
  }
  results$fdr <- p.adjust(results$p.values, "BH")
  results$bon <- p.adjust(results$p.values, "bonferroni")
  cat(paste0(Sys.time(), " End of cor.test.\n"))
  return(results)
}

# Perform correlation analysis
results <- correlation_analysis(data, data2, method, alternative)

# Save results
if (is.null(dataName2)) {
  write.csv(results, paste0(outputDir, "/correlation_for_", dataName1))
} else {
  write.csv(results, paste0(outputDir, "/correlation_for_2_files.csv"))
}

# Function to transform results
cor_tran <- function(data, valueIndex, dataName2) {
  colnames(data)[1] <- "name"
  df <- separate(data = data, col = name, into = c("var1", "var2"), sep = "_vs_")
  unique_vars_1 <- sort(unique(df$var1))
  unique_vars_2 <- sort(unique(df$var2))
  mat <- matrix(NA, ncol = length(unique_vars_1), nrow = length(unique_vars_2),
                dimnames = list(unique_vars_2, unique_vars_1))
  for (i in 1:nrow(df)) {
    row_idx <- which(rownames(mat) == df[i, "var2"])
    col_idx <- which(colnames(mat) == df[i, "var1"])
    mat[row_idx, col_idx] <- df[i, valueIndex]
    if (is.null(dataName2)) {
      mat[col_idx, row_idx] <- df[i, valueIndex]
    }
  }
  return(mat)
}

# Transform and save results
results <- tibble::rownames_to_column(results)

if (!is.null(dataName2)) {
  rvalue_df <- cor_tran(data = results, valueIndex = 3, dataName2 = 2)
  pvalue_df <- cor_tran(data = results, valueIndex = 4, dataName2 = 2)
} else {
  rvalue_df <- cor_tran(data = results, valueIndex = 3, dataName2 = NULL)
  pvalue_df <- cor_tran(data = results, valueIndex = 4, dataName2 = NULL)
}

write.csv(rvalue_df, paste0(outputDir, "/correlation_for_rvalue.csv"))
write.csv(pvalue_df, paste0(outputDir, "/correlation_for_pvalue.csv"))