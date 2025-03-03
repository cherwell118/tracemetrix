#' @title metNor_modified
#' @description Normalize data using different normalization methods.
#' @author Xiaotao Shen,modified by Ziru Chen.
#' @param multiple If multiple = 1, the svr will be built using injection order.
#' If multiple >= 2, the svr is built using top mutiple peaks correlated peaks.
#' For tof data, default is 5, for mrm data, default is 1.
#' @param threads Number of thread.
#' @param path Directory of data
#' @param peakplot Plot feature or not
#' @export
#' @importFrom magrittr %>%
#' @importFrom crayon green
#' @importFrom crayon bgRed
#' @importFrom readr read_csv
#' @importFrom readr cols
#' @import dplyr

metNor_modified <- function(path,threads = 4,multiple=5,peakplot=FALSE){
  ################### MetNormalizer-Data Reading ##################
  cat(crayon::green("Reading data...\n"))
  data <- readr::read_csv(file.path(path, "data.csv"),
                          col_types = readr::cols(),
                          progress = FALSE) %>% as.data.frame()
  
  sample.info <-
    readr::read_csv(file.path(path, "sample.info.csv"), col_types = readr::cols()) %>%
    dplyr::arrange(injection.order)
  
  sample.order <-
    sample.info %>%
    dplyr::filter(class == "Subject") %>%
    dplyr::pull(injection.order) %>%
    as.numeric()
  
  qc.order <-
    sample.info %>%
    dplyr::filter(class == "QC") %>%
    dplyr::pull(injection.order) %>%
    as.numeric()
  
  tags <-
    data %>%
    dplyr::select(-dplyr::one_of(sample.info$sample.name))
  
  tags$name <- gsub("/",".",tags$name)
  
  sample.name <-
    sample.info %>%
    dplyr::filter(class == 'Subject') %>%
    dplyr::pull(sample.name)
  
  qc.name <-
    sample.info %>%
    dplyr::filter(class == 'QC') %>%
    dplyr::pull(sample.name)
  
  sample <-
    data %>%
    dplyr::select(dplyr::one_of(sample.name))
  
  qc <-
    data %>%
    dplyr::select(dplyr::one_of(qc.name))
  
  rownames(sample) <- rownames(qc) <- tags$name
  sample <- t(sample)
  qc <- t(qc)
  tags <- t(tags)
  
  ################### run MetNormalizer ##################
  
  
  cat(crayon::green("SVR normalization...\n"))
  
  MetNormalizer:::SXTsvrNor(
    sample = sample,
    QC = qc,
    tags = tags,
    sample.order = sample.order,
    QC.order = qc.order,
    multiple = multiple,
    path = path,
    rerun = TRUE,
    peakplot = peakplot,
    datastyle = "tof",
    dimension1 = TRUE,
    threads = threads
  )
  cat(crayon::bgRed("All done!\n"))
}

