suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(xcms)))
suppressPackageStartupMessages(suppressWarnings(library(CAMERA)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(pheatmap)))
suppressPackageStartupMessages(suppressWarnings(library(jsonlite)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'sampleInfoName', 'a',1,'character','Name of Sample file.',
  'paraFile','f',1,'character','Json files of parameter list for xcms.',
  'polarity','p',1,'character','polarization mode. positive or negative.'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

if(!is.null(opt$help)){
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
}

# ============================= Main Functions ============================= #

# Function to check and set parameters
check_and_set_parameters <- function(opt) {
  if(is.null(opt$inputDir)) {
    stop("Please set the parameter 'inputDir'")
  }
  inputDir <- opt$inputDir
  
  if(is.null(opt$outputDir)) {
    stop("Please set the parameter 'outputDir'")
  }
  outputDir <- opt$outputDir
  if(!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }
  
  setwd(outputDir)
  
  if(is.null(opt$sampleInfoName)) {
    stop("Please set the parameter 'sampleInfoName'")
  }
  sampleInfoName <- opt$sampleInfoName
  sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
  if(!file.exists(sampleInfoFile)) {
    stop(paste0(sampleInfoFile, " does not exist!"))
  }
  
  if(is.null(opt$paraFile)) {
    stop("The parameter file is missing.")
  }
  paraFile <- opt$paraFile
  paraFilePath <- paste0(inputDir, "/", paraFile)
  paraList <- fromJSON(paraFilePath)
  
  if(!is.null(opt$polarity)) {
    polarity <- opt$polarity
    if(polarity == "positive") {
      adductTable <- NULL
    } else if (polarity == "negative") {
      adductTable <- "negAdducts"
    } else {
      stop("Please check your parameter 'polarity'.")
    }
  }
  
  return(list(inputDir = inputDir, outputDir = outputDir, sampleInfoFile = sampleInfoFile, 
              paraList = paraList, polarity = polarity, adductTable = adductTable))
}

# Function to process sample information
process_sample_info <- function(sampleInfoFile, inputDir) {
  sampleInfo <- read_csv(sampleInfoFile, col_types = cols()) %>% as.data.frame()
  sampleInfo <- sampleInfo[!apply(sampleInfo, 1, function(x) {all(is.na(x))}), ]
  sampleInfo <- sampleInfo[, !apply(sampleInfo, 2, function(x) {all(is.na(x))})]
  if(nrow(sampleInfo) == 0) {
    stop("The sample information file is empty!")
  }
  
  filePaths <- paste0(inputDir, "/", sampleInfo$file.name)
  exist_L <- file.exists(filePaths)
  
  if(!all(exist_L)) {
    if(all(!exist_L)) {
      stop("Can't find any data files, please check your data and sample information files!")
    } else {
      warning("Can't find the following files, please check your data and sample information files!\n",
              paste0(filePaths[!exist_L], "\n"))
      sampleInfo <- sampleInfo[exist_L, ]
      filePaths <- filePaths[exist_L]
    }
  }
  
  return(list(sampleInfo = sampleInfo, filePaths = filePaths))
}

# Function to perform initial data inspection
initial_data_inspection <- function(filePaths, sampleInfo, outputDir) {
  if(anyNA(sampleInfo$classes)) {
    sampleInfo$classes[is.na(sampleInfo$classes)] <- "unknown_group"
  }
  
  pd <- data.frame(class = sampleInfo$classes, sample_name = sampleInfo$sample.names, stringsAsFactors = FALSE)
  raw_data <- readMSData(files = filePaths, pdata = new("NAnnotatedDataFrame", pd), msLevel. = 1, mode = "onDisk")
  
  if(length(raw_data$sample_name) == 0) {
    stop("Your files do not contain any MS1 spectra!")
  } else if(length(raw_data$sample_name) < nrow(sampleInfo)) {
    exist_L <- sampleInfo$sample.names %in% raw_data$sample_name
    warning("The following files do not contain MS1 data.\n", paste0(sampleInfo$file.name[!exist_L], "\n"))
  }
  
  pdfiles <- basename(raw_data@processingData@files)
  sampleInfo <- sampleInfo[sampleInfo$file.name %in% pdfiles, ]
  write.csv(sampleInfo, paste0(outputDir, "/sample.info.csv"))
  
  return(raw_data)
}

# Function to perform chromatographic peak detection
chromatographic_peak_detection <- function(raw_data, paraList, outputDir) {
  para_pd <- paraList$feature_detection
  method <- para_pd$method
  
  if(method == "matchedFilter") {
    para_pd$sigma <- para_pd$FWHM / 2.3548
    paramPeakDetec <- MatchedFilterParam(binSize = para_pd$step, fwhm = para_pd$FWHM, sigma = para_pd$sigma, 
                                         max = para_pd$max_EIC, snthresh = para_pd$snthresh, mzdiff = para_pd$mzdiff)
  } else if(method == "centWave") {
    paramPeakDetec <- CentWaveParam(ppm = para_pd$ppm, peakwidth = c(para_pd$peak_width_min, para_pd$peak_width_max), 
                                    snthresh = para_pd$snthresh, prefilter = c(para_pd$prefilter_peaks, para_pd$prefilter_intensity), 
                                    mzCenterFun = "wMean", integrate = para_pd$integration_method, mzdiff = para_pd$mzdiff, noise = para_pd$noise_filter)
  }
  
  cat(paste0(Sys.time(), " Start: chromatographic peak detection\n"))
  xdata <- findChromPeaks(raw_data, param = paramPeakDetec)
  cat(paste0(Sys.time(), " End: chromatographic peak detection\n"))
  
  para <- paramPeakDetec %>% as.list
  para <- data.frame(Value = unlist(para))
  write.csv(para, paste0(outputDir, "/parameters.csv"))
  
  return(xdata)
}

# Function to perform alignment
alignment <- function(xdata, paraList, outputDir) {
  para_align <- paraList$rtCorrection
  method <- para_align$method
  
  if(method == "obiwarp") {
    paramAdj <- ObiwarpParam(binSize = para_align$profStep)
  } else if(method == "peakgroups") {
    if(para_align$ignore == TRUE) {
      xdata$sample_type <- rep("sg", length(xdata$class))
    } else {
      xdata$sample_type <- xdata$class
    }
    
    pdp <- PeakDensityParam(sampleGroups = xdata$sample_type, bw = para_align$bw, minFraction = 1, minSamples = para_align$minsamp, binSize = para_align$mzwid)
    xdata <- groupChromPeaks(xdata, pdp)
    
    para <- pdp %>% as.list %>% lapply(function(x) if(length(x) > 1) paste(x, collapse = ",") else x)
    para <- data.frame(Value = unlist(para))
    write.table(para, paste0(outputDir, "/parameters.csv"), sep = ",", append = TRUE, col.names = FALSE)
    
    paramAdj <- PeakGroupsParam(minFraction = para_align$minfrac, extraPeaks = para_align$extra, smooth = para_align$smooth, 
                                span = para_align$span, family = para_align$family, subsetAdjust = "average")
  }
  
  cat(paste0(Sys.time(), " Start: adjust retention times\n"))
  xdata <- adjustRtime(xdata, param = paramAdj)
  cat(paste0(Sys.time(), " End: adjust retention times\n"))
  
  para <- paramAdj %>% as.list
  para <- data.frame(Value = unlist(para))
  write.table(para, paste0(outputDir, "/parameters.csv"), sep = ",", append = TRUE, col.names = FALSE)
  
  return(xdata)
}

# Function to perform peak grouping
peak_grouping <- function(xdata, paraList, outputDir) {
  para_gr <- paraList$alignment
  paramGroup <- PeakDensityParam(sampleGroups = xdata$class, bw = para_gr$bw, minFraction = para_gr$minfrac, 
                                 minSamples = para_gr$minsamp, binSize = para_gr$mzwid, maxFeatures = para_gr$max)
  
  para <- paramGroup %>% as.list %>% lapply(function(x) if(length(x) > 1) paste(x, collapse = ",") else x)
  para <- data.frame(Value = unlist(para))
  write.table(para, paste0(outputDir, "/parameters.csv"), sep = ",", append = TRUE, col.names = FALSE)
  
  cat(paste0(Sys.time(), " Start: Peak groups\n"))
  xdata <- groupChromPeaks(xdata, param = paramGroup)
  cat(paste0(Sys.time(), " End: Peak groups\n"))
  
  return(xdata)
}

# Function to perform gap filling
gap_filling <- function(xdata) {
  xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
  return(xdata)
}

# Function to process CAMERA annotation
camera_annotation <- function(xdata, polarity, outputDir) {
  xset <- as(xdata, "xcmsSet")
  xsa <- xsAnnotate(xset)
  xsaF <- groupFWHM(xsa, perfwhm = 0.6)
  xsaC <- groupCorr(xsaF)
  xsaFI <- findIsotopes(xsaC)
  xsaFA <- findAdducts(xsaFI, polarity = polarity)
  peakData <- getPeaklist(xsaFA)
  
  ms1 <- data.frame(name = paste0(rownames(featureData), "_", polarity), mz = Rtmz$mzmed, rt = Rtmz$rtmed, 
                    isotopes = peakData$isotopes, adduct_CAMERA = peakData$adduct)
  write.csv(ms1, paste0(outputDir, "/ms1.peak.table.csv"), row.names = FALSE)
}

# ============================= Main Execution ============================= #

# Step 1: Check and set parameters
params <- check_and_set_parameters(opt)
inputDir <- params$inputDir
outputDir <- params$outputDir
sampleInfoFile <- params$sampleInfoFile
paraList <- params$paraList
polarity <- params$polarity
adductTable <- params$adductTable

# Step 2: Process sample information
sample_info_result <- process_sample_info(sampleInfoFile, inputDir)
sampleInfo <- sample_info_result$sampleInfo
filePaths <- sample_info_result$filePaths

# Step 3: Initial data inspection
raw_data <- initial_data_inspection(filePaths, sampleInfo, outputDir)

# Step 4: Chromatographic peak detection
xdata <- chromatographic_peak_detection(raw_data, paraList, outputDir)

# Step 5: Alignment
xdata <- alignment(xdata, paraList, outputDir)

# Step 6: Peak grouping
xdata <- peak_grouping(xdata, paraList, outputDir)

# Step 7: Gap filling
xdata <- gap_filling(xdata)

# Step 8: CAMERA annotation
camera_annotation(xdata, polarity, outputDir)