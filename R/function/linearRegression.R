suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(car)))
suppressPackageStartupMessages(suppressWarnings(library(visreg)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# Command-line options
commands = matrix(c(
  'help', 'h', 0, "character",'Print this help file and exit.',
  'inputDir', 'i', 1, "character",'Directory of input data.',
  'outputDir', 'o', 1, "character",'Directory for output data.',
  'dataName', 'd', 1, "character",'Information of data file.',
  'y', 'y', 1, "character",'Name of dependent variable.',
  'x', 'x', 1, "character",'Name of independent variable(s).',
  'interaction', 'I', 1, "character",'Interaction between two independent variables.',
  'categorical', 'C', 1, "character",'Names of categorical variable(s), comma separated.'
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
check_required_param("y", opt$y, "Dependent variable name is required.")
check_required_param("x", opt$x, "Independent variable(s) name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
y <- opt$y
x <- opt$x
interaction <- opt$interaction
categorical <- opt$categorical

# Create output directory if it doesn't exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Load data file
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]

# Prepare variables
prepare_variables <- function(data, x, categorical) {
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

vars <- prepare_variables(data, x, categorical)
data <- vars$data
independent <- vars$independent

# Generate formula
generate_formula <- function(y, independent, interaction = NULL) {
  if (!is.null(interaction)) {
    interaction <- unlist(strsplit(interaction, split = ","))
    formula <- as.formula(paste(y, "~", paste(independent, collapse = "+"), "+", paste(interaction, collapse = "+")))
  } else {
    formula <- as.formula(paste(y, "~", paste(independent, collapse = "+")))
  }
  return(formula)
}

fmla <- generate_formula(y, independent, interaction)

# Linear regression analysis
perform_linear_regression <- function(formula, data, outputDir) {
  fit <- lm(formula, data = data)
  
  # Save summary results
  sink_file <- paste0(outputDir, "/summary.txt")
  sink(sink_file)
  print(summary(fit))
  sink()
  
  # Extract and save fit results
  fit_results <- list(coef = coef(summary(fit)), model_stat = summary(fit)$model)
  write.csv(fit_results$coef, paste0(outputDir, "/coefficient_table.csv"), row.names = TRUE)
  write.csv(fit_results$model_stat, paste0(outputDir, "/model_stat_table.csv"), row.names = FALSE)
  
  return(fit)
}

fit <- perform_linear_regression(fmla, data, outputDir)

# Plot regression results
plot_regression_results <- function(fit, data, outputDir, independent) {
  pdf_file <- paste0(outputDir, "/regressionPlot.pdf")
  pdf(pdf_file)
  
  if (length(independent) > 1) {
    visreg::visreg(fit)
  } else if (length(independent) == 1) {
    l <- list(
      a = as.numeric(format(coef(fit)[1], digits = 4)),
      b = as.numeric(format(coef(fit)[2], digits = 4)),
      r2 = format(summary(fit)$r.squared, digits = 4),
      p = format(summary(fit)$coefficients[2, 4], digits = 4)
    )
    eq <- substitute(italic(y) == a + b %.% italic(x)~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
    
    p <- ggplot(data = data, aes(x = data[, independent], y = data[, y])) +
      geom_point() +
      stat_smooth(method = 'lm', formula = y ~ x, colour = 'red') +
      ylab(y) +
      xlab(independent) +
      theme_bw() +
      geom_text(aes(x = (max(data[, independent]) + min(data[, independent])) / 2,
                    y = (max(data[, y]) + min(data[, y])) / 2.5,
                    label = as.character(as.expression(eq))), parse = TRUE)
    print(p)
  }
  dev.off()
}

plot_regression_results(fit, data, outputDir, independent)

# Regression diagnosis plot
pdf(paste0(outputDir, "/diagnosticPlot.pdf"))
par(mfrow = c(2, 2))
plot(fit)
dev.off()

# Outlier detection
pdf(paste0(outputDir, "/influencePlot.pdf"))
car::infIndexPlot(fit)
points <- car::influencePlot(fit)
dev.off()

write.csv(points, paste0(outputDir, "/influencePoint.csv"))
write.csv(car::infIndex(fit), paste0(outputDir, "/influencePoint.all.csv"))