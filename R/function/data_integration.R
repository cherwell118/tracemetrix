suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(batchCorr)))
suppressPackageStartupMessages(suppressWarnings(library(reshape)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'sampleInfoName', 'a',1,"character",'Name of Sample file.',
  'method','m',1,'character','One of median, mean'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

# Check required parameters
if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
}

check_required_param <- function(param_name, param_value, error_message) {
  if (is.null(param_value)) {
    stop(paste0("Please set the parameter '", param_name, "': ", error_message))
  }
}

check_required_param("inputDir", opt$inputDir, "Input directory is required.")
check_required_param("outputDir", opt$outputDir, "Output directory is required.")
check_required_param("dataName", opt$dataName, "Data file name is required.")
check_required_param("sampleInfoName", opt$sampleInfoName, "Sample information file name is required.")
check_required_param("method", opt$method, "Normalization method is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
method <- opt$method

if (!(method %in% c("median", "mean"))) {
  stop("The parameter 'method' must be 'median' or 'mean'.")
}

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Function to load and preprocess data
load_and_preprocess_data <- function(data_file, sample_info_file) {
  data <- read_csv(data_file, col_types = cols(), progress = FALSE) %>% as.data.frame()
  if (nrow(data) == 0) {
    stop("The data file is empty!")
  }
  data <- data[!apply(data, 1, function(x) { all(is.na(x)) }), ]
  data <- data[, !apply(data, 2, function(x) { all(is.na(x)) })]
  rownames(data) <- data$Compound
  compound <- data$Compound

  sample_info <- read_csv(sample_info_file, col_types = cols()) %>% as.data.frame()
  if (nrow(sample_info) == 0) {
    stop("The sample information file is empty!")
  }
  if (is.null(sample_info$batch)) {
    stop("The sample information file has no 'batch' column!")
  }

  return(list(data = data, sample_info = sample_info, compound = compound))
}

# Load data
data_file <- paste0(inputDir, "/", dataName)
sample_info_file <- paste0(inputDir, "/", sampleInfoName)
data_info <- load_and_preprocess_data(data_file, sample_info_file)
data <- data_info$data
sample_info <- data_info$sample_info
compound <- data_info$compound

# Function to prepare data for batch correction
prepare_data_for_batch_correction <- function(data, sample_info) {
  classes <- ifelse("QC" %in% unique(sample_info$classes), c("QC", "Sample"), c("Sample"))
  cal_index <- lapply(classes, function(x) {
    get_index(colnames(data), sample_info, group = x, by = "classes")
  }) %>% unlist() %>% sort()
  index_data <- data.frame(groups = names(cal_index), index = cal_index)
  data_mrs <- data[, -cal_index]
  data <- data[, index_data$index]

  if (length(unique(sample_info$classes)) == 2) {
    data_qc <- data[, sample_info[sample_info$classes == "QC", ]$sample.names]
    rownames(data_qc) <- data$Compound
    data_qc <- tibble::rownames_to_column(data_qc, "Compound")
    data_sample <- data[, sample_info[sample_info$classes == "Sample", ]$sample.names]
    rownames(data_sample) <- data$Compound
    data_sample <- tibble::rownames_to_column(data_sample, "Compound")
  } else {
    data_sample <- data
  }

  return(list(data_mrs = data_mrs, data_qc = data_qc, data_sample = data_sample))
}

# Prepare data
data_prepared <- prepare_data_for_batch_correction(data, sample_info)
data_mrs <- data_prepared$data_mrs
data_qc <- data_prepared$data_qc
data_sample <- data_prepared$data_sample

# Function to plot boxplots before and after batch correction
plot_boxplots <- function(data, sample_info, prefix, output_dir) {
  data_m <- reshape2::melt(data, id = 1)
  batch_num <- unique(sample_info$batch)
  batch_len <- length(batch_num)
  data_m$batch <- "a"

  for (i in 1:batch_len) {
    batch_num_1 <- sample_info$sample.names[which(sample_info$batch %in% batch_num[i])]
    data_m$batch[which(data_m$variable %in% batch_num_1)] <- batch_num[i]
  }

  write.csv(data_m, paste0(output_dir, "/", prefix, "_data.csv"), row.names = FALSE)

  pdf_file <- paste0(output_dir, "/", prefix, "_boxplot.pdf")
  pdf(pdf_file, height = 6, width = 8.5)
  ggplot(data_m, aes(x = variable, y = log10(value), fill = batch)) +
    geom_boxplot(outlier.colour = NA) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15)) +
    xlab("Sample") +
    ylab("log10(Intensity)")
  dev.off()
}

# Plot boxplots before batch correction
plot_boxplots(data_sample, sample_info, "sample_before", outputDir)
if (!is.null(data_qc)) {
  plot_boxplots(data_qc, sample_info[sample_info$classes == "QC", ], "qc_before", outputDir)
}

# Batch correction
data_bf_nor <- data
rownames(data_bf_nor) <- compound
norm_data <- normalizeBatches(peakTableCorr = t(data_bf_nor), 
                              batches = sample_info$batch, 
                              sampleGroup = sample_info$classes, 
                              refGroup = 'QC', 
                              population = 'all')

data_nor <- as.data.frame(t(norm_data$peakTable))
data_nor_out <- cbind(data_mrs, data_nor)
write.csv(data_nor_out, paste0(outputDir, "/data.csv"), row.names = FALSE)
write.csv(sample_info, paste0(outputDir, "/sample_info.csv"), row.names = FALSE)
data_nor <- tibble::rownames_to_column(data_nor, "Compound")

# Plot boxplots after batch correction
plot_boxplots(data_nor, sample_info, "sample_after", outputDir)
if (!is.null(data_qc)) {
  data_nor_qc <- data_nor[, sample_info[sample_info$classes == "QC", ]$sample.names]
  plot_boxplots(data_nor_qc, sample_info[sample_info$classes == "QC", ], "qc_after", outputDir)
}