suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(crayon)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(MetNormalizer)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'dataName', 'd', 1, 'character','Information of data file.',
  'sampleInfoName', 'a',1,'character','Name of Sample file.',
  'method','m',1,'character','Method for missing value process. Including "QC_SVR","median","sum"',
  'peakplot','p',1,'logical','Plot peak or not. Choice: TRUE, FALSE'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

# Help and parameter checks
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

if (!(method %in% c("QC_SVR", "median", "sum"))) {
  stop("The parameter 'method' must be one of 'QC_SVR', 'median', 'sum'.")
}

if (method == "QC_SVR") {
  peakplot <- ifelse(is.null(opt$peakplot), FALSE, as.logical(opt$peakplot))
}

# Create output directory if it doesn't exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Copy sample info file to output directory
sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
if (file.exists(paste0(outputDir, "/", sampleInfoName))) {
  file.remove(paste0(outputDir, "/", sampleInfoName))
}
file.copy(from = sampleInfoFile, to = outputDir)

# Load and preprocess data
load_and_preprocess_data <- function(data_file, sample_info_file) {
  data <- read_csv(data_file, col_types = cols(), progress = FALSE) %>% as.data.frame()
  data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
  data <- data[, !apply(data, 2, function(x) all(is.na(x)))]
  
  sample_info <- read_csv(sample_info_file, col_types = cols()) %>% as.data.frame()
  sample_info <- sample_info[!apply(sample_info, 1, function(x) all(is.na(x))), ]
  sample_info <- sample_info[, !apply(sample_info, 2, function(x) all(is.na(x)))]
  
  return(list(data = data, sample_info = sample_info))
}

data_info <- load_and_preprocess_data(paste0(inputDir, "/", dataName), sampleInfoFile)
data <- data_info$data
sample_info <- data_info$sample_info

# Normalize data
normalize_data <- function(data, sample_info, method, output_dir, peakplot = FALSE) {
  classes <- ifelse("QC" %in% unique(sample_info$classes), c("QC", "Sample"), c("Sample"))
  cal_index <- lapply(classes, function(x) get_index(colnames(data), sample_info, group = x, by = "classes")) %>% unlist()
  
  if (method == "QC_SVR") {
    if (!all(c("injection.order") %in% colnames(sample_info))) {
      stop("Please make sure the column 'injection.order' is in your sample information file.")
    }
    if (!all(c("compound.id", "mz", "rt") %in% colnames(data))) {
      stop("Please make sure the columns 'compound.id', 'mz', and 'rt' are in your data file.")
    }
    
    metnor_dir <- paste0(output_dir, "/QC_SVR_input_data")
    if (!dir.exists(metnor_dir)) {
      dir.create(metnor_dir, recursive = TRUE)
    }
    
    new_data <- cbind(name = data$compound.id, mz = data$mz, rt = data$rt, data[, sort(cal_index)])
    sample_index <- get_index(colnames(new_data), sample_info, group = "Sample", by = "classes")
    sample_data <- new_data[, sample_index]
    sample_sd <- apply(sample_data, 1, sd)
    QC_index <- get_index(colnames(new_data), sample_info, group = "QC", by = "classes")
    QC_data <- new_data[, QC_index]
    QC_sd <- apply(QC_data, 1, sd)
    sd0 <- union(which(QC_sd == 0), which(sample_sd == 0))
    new_data <- new_data[-sd0, ]
    
    write.csv(new_data, paste0(metnor_dir, "/data.csv"), row.names = FALSE)
    
    index <- which(sample_info$classes %in% classes)
    info <- sample_info[index, ]
    info$classes[info$classes == "Sample"] <- "Subject"
    info <- info[c("sample.names", "injection.order", "classes")]
    colnames(info) <- c("sample.name", "injection.order", "class")
    write.csv(info, paste0(metnor_dir, "/sample.info.csv"), row.names = FALSE)
    
    metNor_modified(metnor_dir, threads = 1, multiple = 5, peakplot = peakplot)
    data_filted <- read_csv(paste0(metnor_dir, "/svr_normalization_result/data_svr_normalization.csv"), col_types = cols(), progress = FALSE) %>% as.data.frame()
    data_filted <- data_filted[, -1]
    colnames(data_filted)[c(1, 4, 5)] <- c("compound.id", "sample.rsd", "QC.rsd")
    
    if ("compound.name" %in% colnames(data)) {
      data_new <- data[make.names(data$compound.id) %in% make.names(data_filted$compound.id), ]
      data_filted_new <- data_filted[make.names(data_filted$compound.id) %in% make.names(data$compound.id), ]
      compound_name <- data_new$compound.name
      names(compound_name) <- make.names(data_new$compound.id)
      compound_name <- compound_name[make.names(data_filted$compound.id)]
      data_filted <- cbind(compound.id = data_filted$compound.id, compound.name = compound_name, data_filted[, c(2:ncol(data_filted))])
    }
  } else {
    switch(method,
           median = {
             nor_data <- apply(data[, cal_index], 2, function(x) x / median(x, na.rm = TRUE))
           },
           sum = {
             nor_data <- apply(data[, cal_index], 2, function(x) x / sum(x, na.rm = TRUE))
           })
    data[, cal_index] <- nor_data
    data_filted <- data
    sample_index <- get_index(colnames(data_filted), sample_info, group = "Sample", by = "classes")
    sample <- data_filted[, sample_index]
    data_filted$sample.rsd <- apply(sample, 1, rsdFun)
    
    if ("QC" %in% unique(sample_info$classes)) {
      qc_index <- get_index(colnames(data_filted), sample_info, group = "QC", by = "classes")
      qc <- data_filted[, qc_index]
      data_filted$QC.rsd <- apply(qc, 1, rsdFun)
    }
  }
  return(data_filted)
}

data_filted <- normalize_data(data, sample_info, method, outputDir, peakplot)

# Save normalized data
write.csv(data_filted, paste0(outputDir, "/data.csv"), row.names = FALSE)

# Generate plots
plotDir <- paste0(outputDir, "/plotData")
suppressWarnings(qcplot(inputDir = outputDir, outputDir = plotDir, dataName