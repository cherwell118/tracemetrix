suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(plyr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(globaltest)))
suppressPackageStartupMessages(suppressWarnings(library(qs)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'method', 'm', 1, "character",'Method type in enrichment analysis (ora, qea).',
  'species','s',1,"character","Species of database.",
  'databaseDir','F',1,'character','The directory path of database.',
  'sampleInfoName', 'a', 1, "character",'File containing group information.',
  'databaseName','D',1,"character",'Database name (kegg, smpdb, etc.).',
  'calDataName', 'c', 1, "character",'Data file containing metabolic quantitative data.',
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
check_required_param("method", opt$method, "Method (ora or qea) is required.")
check_required_param("databaseDir", opt$databaseDir, "Database directory is required.")
check_required_param("databaseName", opt$databaseName, "Database name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
method <- opt$method
databaseDir <- opt$databaseDir
databaseName <- opt$databaseName

if (!(method %in% c("ora", "qea"))) {
  stop("The parameter 'method' must be 'ora' or 'qea'.")
}

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Load data files
load_data <- function(file_path, file_name, required_columns = NULL) {
  file <- paste0(file_path, "/", file_name)
  if (!file.exists(file)) {
    stop(paste0(file, " does not exist!"))
  }
  data <- read_csv(file, col_types = cols(), progress = FALSE) %>% as.data.frame()
  data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
  data <- data[, !apply(data, 2, function(x) all(is.na(x)))]
  
  if (!is.null(required_columns)) {
    missing_cols <- required_columns[!required_columns %in% colnames(data)]
    if (length(missing_cols) > 0) {
      stop(paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")))
    }
  }
  return(data)
}

# Load main data
data <- load_data(inputDir, dataName)

if (method == "qea") {
  sampleInfo <- load_data(inputDir, opt$sampleInfoName, required_columns = c("classes"))
  cal <- load_data(inputDir, opt$calDataName, required_columns = c("name"))
}

# Load database
databaseFile <- paste0(databaseDir, "/", opt$species, "_", databaseName, ".rds")
if (!file.exists(databaseFile)) {
  stop(paste0(databaseFile, " does not exist!"))
}
database <- readRDS(databaseFile)

# Enrichment analysis
perform_enrichment_analysis <- function(data, database, method, cal = NULL, sampleInfo = NULL) {
  if (method == "ora") {
    IDlist <- unique(data$kegg_id)
    ID.hits <- database$kegg_id[database$kegg_id %in% IDlist]
    database_hits <- database[database$kegg_id %in% ID.hits, ]
    uniq.count <- length(unique(database$kegg_id))
    
    query_path <- database_hits %>%
      group_by(pathway_names) %>%
      summarise(
        pathway_id = first(pathway_id),
        num = n(),
        compound_names = paste(compound_names, collapse = "; "),
        compound_id = paste(kegg_id, collapse = "; ")
      )
    
    bg_path <- database %>%
      group_by(pathway_names) %>%
      summarise(
        pathway_id = first(pathway_id),
        num = n(),
        compound_names = paste(compound_names, collapse = "; "),
        compound_id = paste(kegg_id, collapse = "; ")
      )
    
    ora_data <- merge(query_path, bg_path, by = "pathway_names")
    ora_data$uniq.count <- uniq.count
    ora_data$q.size <- length(IDlist)
    ora_data$expects <- ora_data$q.size * (ora_data$num / ora_data$uniq.count)
    ora_data$FE <- ora_data$hits_num / ora_data$expects
    ora_data$pvalue <- phyper(ora_data$hits_num - 1, ora_data$num, ora_data$uniq.count - ora_data$num, ora_data$q.size, lower.tail = FALSE)
  } else if (method == "qea") {
    hits.list <- split(database$compound_names, database$pathway_names)
    classes <- as.numeric(factor(sampleInfo$classes))
    gt <- globaltest::gt(classes, cal, subsets = hits.list)
    results <- globaltest::result(gt)
    ora_data <- data.frame(
      pathway = names(hits.list),
      hits = sapply(hits.list, length),
      set.num = sapply(hits.list, function(x) length(unique(x))),
      Statistic_Q = results$statistic,
      Expectes_Q = results$expected,
      pvalue = results$p.value
    )
  }
  
  ora_data$holm <- p.adjust(ora_data$pvalue, "holm")
  ora_data$fdr <- p.adjust(ora_data$pvalue, "fdr")
  ora_data$bonf <- p.adjust(ora_data$pvalue, "bonf")
  ora_data <- ora_data[order(ora_data$pvalue), ]
  return(ora_data)
}

# Perform enrichment analysis
ora_data <- perform_enrichment_analysis(data, database, method, cal, sampleInfo)

# Save results
write.csv(ora_data, paste0(outputDir, "/enrich_data.csv"), row.names = FALSE)

# Plotting
plot_enrichment_results <- function(ora_data, output_dir, method, n = 20) {
  plot_data <- ora_data[1:n, ]
  
  if (method == "ora") {
    p1 <- ggplot(plot_data, aes(x = FE, y = pathway_names)) +
      geom_point(aes(size = hits_num, color = -log10(pvalue))) +
      scale_color_gradient(low = "green", high = "red") +
      labs(color = expression(-log[10](pvalue)), size = "Hits", title = "Pathway Enrichment Dotplot") +
      ylab("") + xlab("Fold Enrichment") + theme_bw()
    
    p2 <- ggplot(plot_data, aes(x = pathway_names, y = FE, fill = pvalue)) +
      geom_bar(stat = 'identity', color = 'black', width = 0.65) +
      coord_flip() +
      scale_fill_gradient(low = 'red', high = 'darkgoldenrod1') +
      labs(title = "Pathway Enrichment Histogram") +
      ylab("Fold Enrichment") + xlab("Pathway Name") + theme_bw()
  } else if (method == "qea") {
    p1 <- ggplot(plot_data, aes(x = Statistic_Q, y = pathway_names)) +
      geom_point(aes(size = hits_num, color = -log10(pvalue)))