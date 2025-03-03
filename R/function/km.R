suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(survival)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(survminer)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'time','t',1, "character",'Survival time of patients.',
  'event','e',1, "character",'Survival event of patients.',
  'groups','g',1, "character",'Groups of patients.'
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
check_required_param("time", opt$time, "Survival time column name is required.")
check_required_param("event", opt$event, "Survival event column name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
time <- opt$time
event <- opt$event
groups <- opt$groups

# Create output directory if it doesn't exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Load data file
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]

# Perform Kaplan-Meier analysis
perform_km_analysis <- function(data, time, event, groups = NULL) {
  if (!is.null(groups)) {
    fmla <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", groups))
  } else {
    fmla <- as.formula(paste0("Surv(", time, ",", event, ") ~ 1"))
  }
  
  fit <- survfit(fmla, data = data)
  fit$call$formula <- fmla
  
  return(fit)
}

fit <- perform_km_analysis(data, time, event, groups)

# Save summary results
sink_file <- paste0(outputDir, "/km.summary.txt")
sink(sink_file)
print(fit)
cat("\n")
print(summary(fit))
sink()

# Extract and save median survival time and survival table
extract_and_save_results <- function(fit, outputDir, groups = NULL) {
  summary <- summary(fit)
  
  # Median survival time
  mst <- summary$table
  if (!is.matrix(mst)) {
    write.csv(t(mst), paste0(outputDir, "/km.mst.csv"), row.names = FALSE)
  } else {
    write.csv(mst, paste0(outputDir, "/km.mst.csv"))
  }
  
  # Survival table
  if (!is.null(groups)) {
    surv_table <- data.frame(
      strata = summary$strata,
      time = summary$time,
      risk = summary$n.risk,
      event = summary$n.event,
      survival = summary$surv,
      std.err = summary$std.err,
      lower95CI = summary$lower,
      upper95CI = summary$upper
    )
  } else {
    surv_table <- data.frame(
      time = summary$time,
      risk = summary$n.risk,
      event = summary$n.event,
      survival = summary$surv,
      std.err = summary$std.err,
      lower95CI = summary$lower,
      upper95CI = summary$upper
    )
  }
  write.csv(surv_table, paste0(outputDir, "/km.survival.table.csv"), row.names = FALSE)
}

extract_and_save_results(fit, outputDir, groups)

# Plot Kaplan-Meier curves
plot_km_curves <- function(fit, data, outputDir, groups = NULL) {
  pdf_file <- paste0(outputDir, "/km.Plot.pdf")
  pdf(pdf_file)
  
  if (!is.null(groups)) {
    ggsurvplot(
      fit,
      data = data,
      pval = TRUE,
      conf.int = TRUE,
      linetype = "strata",
      surv.median.line = "hv",
      ggtheme = theme_bw(),
      title = "Kaplan-Meier Curves"
    )
  } else {
    ggsurvplot(
      fit,
      data = data,
      conf.int = TRUE,
      surv.median.line = "hv",
      ggtheme = theme_bw(),
      title = "Kaplan-Meier Curves"
    )
  }
  dev.off()
}

plot_km_curves(fit, data, outputDir, groups)

# Perform survival difference test
if (!is.null(groups)) {
  perform_survdiff <- function(fit, data, outputDir, groups) {
    if (nlevels(factor(data[, groups])) == 2) {
      surv_diff <- survdiff(fit$call$formula, data = data)
      n <- surv_diff$n
      obs <- surv_diff$obs
      exp <- surv_diff$exp
      oee <- ((obs - exp)^2) / exp
      oev <- (((obs - exp)^2) / diag(surv_diff$var))
      df <- (sum(1 * (exp > 0))) - 1
      chisq <- surv_diff$chisq
      p <- format.pval(pchisq(chisq, df, lower.tail = FALSE))
      
      diff_table <- cbind(n, obs, exp, oee, oev)
      dimnames(diff_table) <- list(names(surv_diff$n), c("N", "Observed", "Expected", "(O-E)^2/E", "(O-E)^2/V"))
      write.csv(diff_table, paste0(outputDir, "/statistics.csv"))
      
      cat(paste0("Chisq=", format(round(chisq, 1)), " on", df, "degrees of freedom, p=", p),
          file = paste0(outputDir, "/statistics.txt"))
    }
  }
  
  perform_survdiff(fit, data, outputDir, groups)
}