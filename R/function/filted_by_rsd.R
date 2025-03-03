suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(pheatmap)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "logical",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'sampleInfoName','a',1,"character",'Name of Sample file.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'rsd','q',1,'numeric','the rate of qc.rsd'
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

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName

# Set default RSD value if not provided
rsd <- ifelse(is.null(opt$rsd), 30, as.numeric(opt$rsd))

# Create output directory if it doesn't exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Copy sample info file to output directory
sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
outputSampleInfoFile <- paste0(outputDir, "/", sampleInfoName)

if (file.exists(outputSampleInfoFile)) {
  file.remove(outputSampleInfoFile)
}
file.copy(from = sampleInfoFile, to = outputSampleInfoFile)

# Load data files
load_data <- function(file_path, file_name) {
  file <- paste0(file_path, "/", file_name)
  if (!file.exists(file)) {
    stop(paste0(file, " does not exist!"))
  }
  data <- read_csv(file, col_types = cols(), progress = FALSE) %>% as.data.frame()
  data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
  data <- data[, !apply(data, 2, function(x) all(is.na(x)))]
  return(data)
}

data <- load_data(inputDir, dataName)
sampleInfo <- load_data(inputDir, sampleInfoName)

# Filter data based on RSD
filter_data_by_rsd <- function(data, rsd_threshold) {
  if ("QC" %in% unique(sampleInfo$classes)) {
    nrow_before <- nrow(data)
    cat(paste0("The feature size before filtering is ", nrow_before, ".\n"))
    data_filted <- data[data$QC.rsd <= rsd_threshold, ]
    nrow_after <- nrow(data_filted)
    cat(paste0("The feature size after filtering is ", nrow_after, ".\n"))
    cat(paste0(nrow_before - nrow_after, " features have been removed.\n"))
    
    if (nrow(data_filted) == 0) {
      stop("Data after filtering is empty.")
    }
    return(data_filted)
  } else {
    warning("Only QC.rsd can be used for filtering data.")
    return(data)
  }
}

data_filted <- filter_data_by_rsd(data, rsd)

# Save filtered data
write.csv(data_filted, paste0(outputDir, "/data.csv"), row.names = FALSE)
write.csv(sampleInfo, paste0(outputDir, "/sample_info.csv"), row.names = FALSE)

# Generate QC plots
plotDir <- paste0(outputDir, "/plotData")
if (!dir.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

suppressWarnings(qcplot(inputDir = outputDir, outputDir = plotDir, dataName = "data.csv", sampleInfoName = "sample_info.csv"))