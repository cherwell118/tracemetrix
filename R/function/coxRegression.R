suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(survival)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(survminer)))
suppressPackageStartupMessages(suppressWarnings(library(rms)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(car)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'time','t',1, "character",'Survival time of patients.',
  'event','e',1, "character",'Survival event of patients.',
  'independent','I',1, "character",'Name of independent variant(s).',
  'unorderC','U',1,"character",'Names of unorder categorical variable(s),comma separated.',
  'orderC','D',1,"character",'Names of order categorical variable(s),comma separated.'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
}

# Extract options
inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
time <- opt$time
event <- opt$event
independent <- opt$independent
unorderC <- opt$unorderC
orderC <- opt$orderC

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Load data
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) { all(is.na(x)) }), ]
data <- data[, !apply(data, 2, function(x) { all(is.na(x)) })]

# Function to preprocess categorical variables
preprocess_categorical <- function(data, unorderC = NULL, orderC = NULL) {
  if (!is.null(unorderC)) {
    unorderC <- unlist(strsplit(unorderC, split = ","))
    for (variable in unorderC) {
      data[[variable]] <- as.factor(data[[variable]])
    }
  }
  
  if (!is.null(orderC)) {
    orderC <- unlist(strsplit(orderC, split = ","))
    for (variable in orderC) {
      data[[variable]] <- factor(data[[variable]], ordered = TRUE)
    }
  }
  return(data)
}

# Preprocess data
data <- preprocess_categorical(data, unorderC, orderC)

# Function to generate Cox regression formula
generate_formula <- function(time, event, independent) {
  independent <- unlist(strsplit(independent, split = ","))
  formula <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", paste(independent, collapse = "+")))
  return(formula)
}

# Generate formula
formula <- generate_formula(time, event, independent)

# Function to perform Cox regression and save results
cox_regression_analysis <- function(data, formula, outputDir) {
  fit <- coxph(formula, data = data)
  fit$call$formula <- formula
  
  sink(paste0(outputDir, "/coxph.summary.txt"))
  print(summary(fit))
  sink()
  
  summary <- summary(fit)
  coef.int <- summary$conf.int
  colnames(coef.int) <- c("exp(coef)", "exp(-coef)", paste0("exp(coef)_", c("lower .95", "upper .95")))
  
  summTable <- cbind(coef = summary$coefficients[, 1], coef.int, summary$coefficients[, -c(1, 2)])
  write.csv(summTable, paste0(outputDir, "/conf.statistics.csv"))
  
  sink(paste0(outputDir, "/coxph.statistics.txt"))
  if (!is.null(summary$concordance)) {
    cat("Concordance=", format(round(summary$concordance[1], 3)), " (se =", format(round(summary$concordance[2], 3)), ")\n")
  }
  cat("Likelihood ratio test= ", format(round(summary$logtest["test"], 2)), "  on ", summary$logtest["df"], " df,", "   p=", format.pval(summary$logtest["pvalue"], digits = 3), "\n", sep = "")
  cat("Wald test            = ", format(round(summary$waldtest["test"], 2)), "  on ", summary$waldtest["df"], " df,", "   p=", format.pval(summary$waldtest["pvalue"], digits = 3), "\n", sep = "")
  cat("Score (logrank) test = ", format(round(summary$sctest["test"], 2)), "  on ", summary$sctest["df"], " df,", "   p=", format.pval(summary$sctest["pvalue"], digits = 3), sep = "")
  if (!is.null(summary$robscore)) {
    cat(",   Robust = ", format(round(summary$robscore["test"], 2)), "  p=", format.pval(summary$robscore["pvalue"], digits = 3), "\n\n", sep = "")
  }
  sink()
}

# Perform Cox regression analysis
cox_regression_analysis(data, formula, outputDir)

# Function to plot ZPH test results
plot_zph <- function(fit, outputDir) {
  zph <- cox.zph(fit)$table
  write.csv(zph, paste0(outputDir, "/czph.csv"))
  
  pdf(paste0(outputDir, "/zph.Plot.pdf"))
  plot(cox.zph(fit))
  dev.off()
}

# Plot ZPH test results
fit <- coxph(formula, data = data)
plot_zph(fit, outputDir)

# Function to create Forest plot
create_forest_plot <- function(fit, data, outputDir) {
  n <- nrow(coef(fit)) + 1
  pdf(paste0(outputDir, "/forest.Plot.pdf"), height = n * 1.3, width = n * 1.3 * 1.5)
  ggforest(fit, data = data, main = "Hazard Ratio", cpositions = c(0.02, 0.22, 0.4), fontsize = 0.8, refLabel = "reference", noDigits = 2)
  dev.off()
}

# Create Forest plot
create_forest_plot(fit, data, outputDir)

# Function to create Nomogram
create_nomogram <- function(data, formula, outputDir) {
  ddata <- datadist(data)
  options(datadist = "ddata")
  
  fit2 <- cph(formula, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36, data = data)
  fit2$call$formula <- formula
  
  surv <- Survival(fit2)
  nom <- nomogram(fit2, fun = list(function(x) surv(36, x), function(x) surv(60, x), function(x) surv(120, x)), 
                  lp = FALSE, funlabel = c("3-year survival", "5-year survival", "10-year survival"), maxscale = 100, 
                  fun.at = list(c(0.95, 0.94), c(0.97, 0.95, 0.91, 0.9, 0.88), c(0.95, 0.9, 0.85, 0.8, 0.75)))
  
  pdf(paste0(outputDir, "/nomogram.Plot.pdf"))
  plot(nom)
  dev.off()
}

# Create Nomogram
create_nomogram(data, formula, outputDir)