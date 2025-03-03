# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(metID)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'polarity', 'p', 1, 'character','Polarization mode: positive or negative.',
  'annotationType', 't', 1, 'character','Annotation type: ms1, ms2 or all.',
  'databaseDir', 'F', 1, 'character','The directory path of database.',
  'database', 'd', 1, 'character','Database name(s) chosen by user, comma separated.',
  'ms1.match.ppm', 's', 1, 'integer','MS1 tolerance for database match, ppm.',
  'candidateNum', 'c', 1, 'integer','Number of candidate annotations.',
  'column', 'u', 1, 'character','Column type: hilic or rp.',
  'ms1file', 'm', 1, 'character','Name of MS1 file.',
  'ms2files', 'f', 1, 'character','Names of MS2 file(s), comma separated.'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

# 帮助信息
if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
}

# 参数检查
check_required_param <- function(param_name, param_value, error_message) {
  if (is.null(param_value)) {
    stop(paste0("Please set the parameter '", param_name, "': ", error_message))
  }
}

check_required_param("inputDir", opt$inputDir, "Input directory is required.")
check_required_param("outputDir", opt$outputDir, "Output directory is required.")
check_required_param("polarity", opt$polarity, "Polarization mode is required.")
check_required_param("database", opt$database, "Database name(s) are required.")
check_required_param("databaseDir", opt$databaseDir, "Database directory is required.")
check_required_param("ms1file", opt$ms1file, "MS1 file name is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
polarity <- opt$polarity
database <- unlist(strsplit(opt$database, split = ","))
databaseDir <- opt$databaseDir
ms1file <- opt$ms1file
annotationType <- ifelse(is.null(opt$annotationType), "all", opt$annotationType)
ms1.match.ppm <- ifelse(is.null(opt$ms1.match.ppm), 15, as.numeric(opt$ms1.match.ppm))
candidateNum <- ifelse(is.null(opt$candidateNum), 3, as.numeric(opt$candidateNum))
column <- ifelse(is.null(opt$column), "rp", opt$column)

# 创建输出目录
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# 复制文件到输出目录
copy_files_to_output <- function(file_name, inputDir, outputDir) {
  if ((!file.exists(paste0(outputDir, "/", file_name))) & file.exists(paste0(inputDir, "/", file_name))) {
    file.copy(from = paste0(inputDir, "/", file_name), to = outputDir, overwrite = TRUE, recursive = TRUE)
  } else if ((!file.exists(paste0(outputDir, "/", file_name))) & !(file.exists(paste0(inputDir, "/", file_name)))) {
    stop(paste0("Cannot find file ", file_name, "!"))
  }
}

copy_files_to_output(ms1file, inputDir, outputDir)

if ((annotationType == "ms2") | (annotationType == "all")) {
  if (!is.null(opt$ms2files)) {
    ms2files <- unlist(strsplit(opt$ms2files, split = ","))
    for (ms2file in ms2files) {
      copy_files_to_output(ms2file, inputDir, outputDir)
    }
  } else {
    stop("MS2 file is missing, please check your parameter 'ms2files'.")
  }
}

# 准备数据库
prepare_database <- function(databaseDir, database, annotationType) {
  database_info <- read_csv(paste0(databaseDir, "/db_info.csv"), col_types = cols(), progress = FALSE) %>% as.data.frame()
  rownames(database_info) <- database_info$db_name
  
  if (annotationType == "ms1") {
    databases <- (database_info[database, ] %>% filter(type == "MS1"))$full_name
    if (length(databases) == 0) {
      stop("No MS1 database is selected. Please check your parameter 'database'.")
    }
  } else if (annotationType == "ms2") {
    databases <- (database_info[database, ] %>% filter(type == "MS2"))$full_name
    if (length(databases) == 0) {
      stop("No MS2 database is selected. Please check your parameter 'database'.")
    }
  } else if (annotationType == "all") {
    databases <- database_info[database, ]$full_name
    if (length(databases) == 0) {
      stop("Neither MS1 nor MS2 database is selected. Please check your parameter 'database'.")
    }
  }
  
  return(databases)
}

databases <- prepare_database(databaseDir, database, annotationType)

# 准备MS1文件
ms1_data <- read_csv(paste0(inputDir, "/", ms1file), col_types = cols(), progress = FALSE) %>% as.data.frame()
ms1_file <- data.frame(name = ms1_data$compound.id, mz = ms1_data$mz, rt = ms1_data$rt)
write_csv(ms1_file, paste0(outputDir, "/ms1_file.csv"))

# MS1注释
annotate_ms1 <- function(databases, ms1.match.ppm, polarity, column, candidateNum, outputDir) {
  annotate_result <- lapply(databases, function(database) {
    cat(paste0(Sys.time(), " Start running ", database, "\n"))
    annotate_result_temp <- identify_metabolites(
      ms1.data = "ms1_file.csv",
      ms1.match.ppm = ms1.match.ppm,
      rt.match.tol = 1000000,
      polarity = polarity,
      column = column,
      path = outputDir,
      total.score.tol = 0.2,
      candidate.num = candidateNum,
      database = database,
      threads = 5
    )
    cat(paste0(Sys.time(), " End of running ", database, "\n"))
    annotation_table <- get_identification_table_m(annotate_result_temp, candidate.num = candidateNum, type = "new")
    annotation_table %<>% filter(!is.na(Compound.name)) %>% mutate(Level = 3) %>% as.data.frame
    return(annotation_table)
  })
  
  names(annotate_result) <- databases
  annotate_result_ms1 <- do.call(rbind, annotate_result)
  return(annotate_result_ms1)
}

if ((annotationType == "ms1") | (annotationType == "all")) {
  ms1_databases <- (database_info[database, ] %>% filter(type == "MS1"))$full_name
  if (length(ms1_databases) > 0) {
    file.copy(from = paste0(databaseDir, "/", ms1_databases), to = outputDir, overwrite = TRUE, recursive = TRUE)
    annotate_result_ms1 <- annotate_ms1(ms1_databases, ms1.match.ppm, polarity, column, candidateNum, outputDir)
    file.remove(paste0(outputDir, "/", ms1_databases))
  }
}

# MS2注释
annotate_ms2 <- function(databases, ms1.match.ppm, polarity, column, candidateNum, ms2files, outputDir) {
  annotate_result <- lapply(databases, function(database) {
    cat(paste0(Sys.time(), " Start running ", database, "\n"))
    annotate_result_temp <- identify_metabolites(
      ms1.data = "ms1_file.csv",
      ms2.data = ms2files,
      ms2.match.tol = 0.1,
      ce = "all",
      ms1.match.ppm = ms1.match.ppm,
      rt.match.tol = 1000000,
      polarity = polarity,
      column = column,
      path = outputDir,
      ms2.match.ppm = 50,
      total.score.tol = 0.2,
      candidate.num = candidateNum,
      database = database,
      threads = 5
    )
    cat(paste0(Sys.time(), " End of running ", database, "\n"))
    annotation_table <- get_identification_table_m(annotate_result_temp, candidate.num = candidateNum, type = "new")
    if (!is.null(annotation_table)) {
      annotation_table %<>% filter(!is.na(Compound.name)) %>% mutate(Level = 2) %>% as.data.frame
    }
    return(annotation_table)
  })
  
  names(annotate_result) <- databases
  annotate_result_ms2 <- do.call(rbind, annotate_result)
  return(annotate_result_ms2)
}

if ((annotationType == "ms2") | (annotationType == "all")) {
  ms2_databases <- (database_info[database, ] %>% filter(type == "MS2"))$full_name
  if (length(ms2_databases) > 0) {
    file.copy(from = paste0(databaseDir, "/", ms2_databases), to = outputDir, overwrite = TRUE, recursive = TRUE)
    annotate_result_ms2 <- annotate_ms2(ms2_databases, ms1.match.ppm, polarity, column, candidateNum, ms2files, outputDir)
    file.remove(paste0(outputDir, "/", ms2_databases))
  }
}

# 整合注释结果
integrate_annotations <- function(annotate_result_ms1, annotate_result_ms2, annotationType) {
  if (annotationType == "ms1") {
    annotate_result <- annotate_result_ms1
  } else if (annotationType == "ms2") {
    annotate_result <- annotate_result_ms2
  } else if (annotationType == "all") {
    annotate_result <- rbind(annotate_result_ms1, annotate_result_ms2)
  }
  return(annotate_result)
}

annotate_result <- integrate_annotations(annotate_result_ms1, annotate_result_ms2, annotationType)

# 清理注释结果
clean_annotation <- function(annotate_result) {
  # 删除不必要的列
  cols_to_remove <- c("MS2.spectra.name", "Candidate.number", "SS", "CE", "mz.match.score", "RT.error", "RT.match.score")
  showFile <- annotate_result[, setdiff(names(annotate_result), cols_to_remove)]
  
  # 替换数据库名称
  replacement <- c("bloodExposomeMS1Database_1.0", "drugbankMS1Database5.1.8", "hmdbMS1Database0.0.1",
                  "keggMS1Database_1.0", "metlinMS1Database2021", "t3dbMS1Database_1.0",
                  "hmdbDatabase0.0.2", "massbankDatabase0.0.2", "metlinMS2Database2021",
                  "monaDatabase0.0.2")
  replacement_values <- c("BloodExposome", "DrugBank", "HMDB", "KEGG", "Metlin", "T3DB",
                          "HMDB", "MassBank", "Metlin", "MONA")
  
  showFile$Database <- replacement_values[match(showFile$Database, replacement)]
  
  # 选择最佳注释结果
  showFile <- showFile[order(showFile$name, decreasing = FALSE), ]
  finalFile <- showFile[-c(1:nrow(showFile)), ]
  colnames(finalFile) <- colnames(showFile)
  
  featureIndex <- sort(unique(showFile$name))
  for (i in 1:length(featureIndex)) {
    df <- showFile[which(showFile$name %in% featureIndex[i]), ]
    if ("ms2" %in% unique(df$Level)) {
      df_ms2 <- df[which(df$Level == "ms2"), ]
      if (nrow(df_ms2) > 1) {
        if (any(grepl("HMDB", df_ms2$HMDB.ID))) {
          df_ms2_hmdb <- df_ms2[grepl("HMDB", df_ms2$HMDB.ID), ]
          if (nrow(df_ms2_hmdb) > 1) {
            max_score <- max(df_ms2_hmdb$Total.score)
            df_ms2_hmdb_max <- df_ms2_hmdb[df_ms2_hmdb$Total.score == max_score, ]
            finalFile <- rbind(finalFile, df_ms2_hmdb_max[1, ])
          } else {
            finalFile <- rbind(finalFile, df_ms2_hmdb)
          }
        } else {
          max_score <- max(df_ms2$Total.score)
          df_ms2_max <- df_ms2[df_ms2$Total.score == max_score, ]
          finalFile <- rbind(finalFile, df_ms2_max[1, ])
        }
      } else {
        finalFile <- rbind(finalFile, df_ms2)
      }
    } else {
      df_ms1 <- df[which(df$Level == "ms1"), ]
      if (any(grepl("HMDB", df_ms1$HMDB.ID))) {
        df_ms1_hmdb <- df_ms1[grepl("HMDB", df_ms1$HMDB.ID), ]
        if (nrow(df_ms1_hmdb) > 1) {
          max_score <- max(df_ms1_hmdb$Total.score)
          df_ms1_hmdb_max <- df_ms1_hmdb[df_ms1_hmdb$Total.score == max_score, ]
          finalFile <- rbind(finalFile, df_ms1_hmdb_max[1, ])
        } else {
          finalFile <- rbind(finalFile, df_ms1_hmdb)
        }
      } else {
        max_score <- max(df_ms1$Total.score)
        df_ms1_max <- df_ms1[df_ms1$Total.score == max_score, ]
        finalFile <- rbind(finalFile, df_ms1_max[1, ])
      }
    }
  }
  
  finalFile$mz.error <- round(as.numeric(finalFile$mz.error), digits = 6)
  finalFile[is.na(finalFile)] <- ""
  return(finalFile)
}

if (!is.null(annotate_result)) {
  write_csv(annotate_result, paste0(outputDir, "/annotation_table_", annotationType, ".csv"))
  
  annotation_table_s <- clean_annotation(annotate_result)
  
  new_cols <- c("name", "mz", "rt", "Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID", "Adduct", "Lab.ID", "mz.error", "Total.score", "Database", "Level")
  annotation_table_s <- annotation_table_s[, new_cols]
  
  write_csv(annotation_table_s, paste0(outputDir, "/annotation_table_", annotationType, "_s.csv"))
  
  annotate_name <- annotate_result %>% select(c(name, Compound.name))
  annotate_name <- split(annotate_name$Compound.name, annotate_name$name) %>%
    lapply(function(x) { paste0(x, collapse = ",") }) %>% data.frame()
  
  annotate_name <- data.frame(compound.id = names(annotate_name), annotate.name = t(annotate_name))
  write_csv(annotate_name, paste0(outputDir, "/annotate_name_", annotationType, ".csv"))
  
  ms1_data <- ms1_data %>% left_join(annotate_name, by = "compound.id")
  ms1_data$compound.name <- ms1_data$annotate.name
  ms1_data <- ms1_data %>% select(-annotate.name)
  write_csv(ms1_data, paste0(outputDir, "/featureData_anno_", annotationType, ".csv"))
} else {
  warning("No feature is annotated!")
}

# 重命名MS2输出目录
ms2FileRename <- function(fileDir, oldName, newName) {
  old <- paste0(fileDir, "/", oldName)
  new <- paste0(fileDir, "/", newName)
  if (file.exists(old)) {
    file.rename(old, new)
  } else {
    cat("No database directory exists\n")
  }
}

ms2FileRename(fileDir = outputDir, oldName = "ms2plot_massbankDatabase0.0.2", newName = "ms2plot_MassBank")
ms2FileRename(fileDir = outputDir, oldName = "ms2plot_metlinMS2Database2021", newName = "ms2plot_Metlin")
ms2FileRename(fileDir = outputDir, oldName = "ms2plot_monaDatabase0.0.2", newName = "ms2plot_MONA")
ms2FileRename(fileDir = outputDir, oldName = "ms2plot_hmdbDatabase0.0.2", newName = "ms2plot_HMDB")

cat(paste0(Sys.time(), " Ending running Annotation\n"))