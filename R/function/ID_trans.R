suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(plyr)))
suppressPackageStartupMessages(suppressWarnings(library(qs)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "logical",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'inputType', 'm', 1, "character",'Input file type: compound list (CList) or quantitative data (QData).',
  'databaseDir','F',1,'character','The directory path of database.',
  'Type', 't', 1, "character",'Type of metabolite ID: kegg, hmdb, names, metlin, chebi, pubchem.',
  'dataName', 'd', 1, "character",'Information of data file.'
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
check_required_param("Type", opt$Type, "Type of metabolite ID is required.")
check_required_param("databaseDir", opt$databaseDir, "Database directory is required.")
check_required_param("inputType", opt$inputType, "Input type (CList or QData) is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
Type <- opt$Type
databaseDir <- opt$databaseDir
inputType <- opt$inputType

# Create output directory if it doesn't exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

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
databaseFile <- paste0(databaseDir, "/compound_db.qs")
database <- qread(databaseFile)
database$name <- iconv(database$name, "WINDOWS-1252", "UTF-8")

# Initialize variables
if (inputType == "CList") {
  cpd.ori <- data$name
} else if (inputType == "QData") {
  cpd.ori <- names(data)[-1]
} else {
  stop("Invalid input type. Use 'CList' for compound list or 'QData' for quantitative data.")
}

# Mapping function
map_metabolites <- function(cpd.ori, database, Type) {
  hit.inx <- match(tolower(cpd.ori), tolower(database[[Type]]))
  match.values <- database$name[hit.inx]
  match.state <- as.numeric(!is.na(hit.inx))
  
  name_map <- database[hit.inx, ]
  name_map <- name_map[!is.na(hit.inx), ]
  
  names.map.out <- cbind(cpd.ori, name_map[, c("name", "hmdb_id", "pubchem_id", "chebi_id", "kegg_id", "metlin_id", "lipid")])
  colnames(names.map.out) <- c("query", "name", "hmdb_id", "pubchem_id", "chebi_id", "kegg_id", "metlin_id", "lipid")
  
  return(list(names.map.out = names.map.out, match.state = match.state))
}

# Perform mapping
mapping_result <- map_metabolites(cpd.ori, database, Type)
names.map.out <- mapping_result$names.map.out
match.state <- mapping_result$match.state

# Save mapped data
write.csv(names.map.out, paste0(outputDir, "/name_map.csv"), row.names = FALSE)

# Handle fuzzy matching for unmatched compounds
if (Type == "names") {
  value.diff <- setdiff(cpd.ori, na.omit(names.map.out$name))
  metabolites.fuzzy <- lapply(value.diff, function(x) {
    agrep(x, database$name, ignore.case = TRUE, value = TRUE, max.distance = 0.05, useBytes = TRUE)
  })
  
  metabolites.fuzzy <- metabolites.fuzzy[sapply(metabolites.fuzzy, length) > 0]
  
  if (length(metabolites.fuzzy) > 0) {
    DiffOutputDir <- paste0(outputDir, "/Diff")
    if (!dir.exists(DiffOutputDir)) {
      dir.create(DiffOutputDir, recursive = TRUE)
    }
    
    for (j in seq_along(metabolites.fuzzy)) {
      hit.diff.inx <- match(tolower(metabolites.fuzzy[[j]]), tolower(database$name))
      diff.values <- database[hit.diff.inx, ]
      diff.values <- cbind(diff.values[, 2], diff.values[, 1], diff.values[, 3:7])
      colnames(diff.values) <- colnames(names.map.out)[-1]
      write.csv(diff.values, paste0(DiffOutputDir, "/", names(metabolites.fuzzy)[j], ".csv"), row.names = FALSE)
    }
  }
}

# Save transposed data if input type is QData
if (inputType == "QData") {
  data <- data[, colnames(data) %in% cpd.ori[match.state == 1]]
  data <- data.frame(t(data))
  data <- tibble::rownames_to_column(data, var = "name")
  write.csv(data, paste0(outputDir, "/trans_data.csv"), row.names = FALSE)
}