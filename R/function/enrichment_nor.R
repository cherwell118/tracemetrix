
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(crayon)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
#suppressPackageStartupMessages(suppressWarnings(library(MetNormalizer)))

commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'dataName', 'd', 1, 'character','Information of data file.',
  'sampleInfoName', 'a',1,'character','Name of Sample file.',
  'method','m',1,'character','Method for missing value process. Including ,"median","sum"'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

if(!is.null(opt$help)){
  cat(paste(getopt(commands, usage = T), "\n"))
  q(status=1)
}

if(is.null(opt$inputDir)){
  stop("Please set the parameter 'inputDir'")
}else{
  inputDir <- opt$inputDir
}

if(is.null(opt$outputDir)){
  stop("Please set the parameter 'outputDir'")
}else{
  outputDir <- opt$outputDir
  if(!dir.exists(outputDir)){
    dir.create(outputDir,recursive = TRUE)
  }
}

if(is.null(opt$dataName)){
  stop("Please set the parameter 'dataName'")
}else{
  dataName <- opt$dataName
  dataFile <- paste0(inputDir,"/",dataName)
  if(!file.exists(dataFile)){
    stop(paste0(dataFile," do not exist!"))
  }
}
if(is.null(opt$sampleInfoName)){
  stop("Please set the parameter 'sampleInfoName'")
}else{
  sampleInfoName <- opt$sampleInfoName
  sampleInfoFile <- paste0(inputDir,"/",sampleInfoName)
  if(!file.exists(sampleInfoFile)){
    stop(paste0(sampleInfoFile," do not exist!"))
  }
}
if(is.null(opt$method)){
  stop("Please set the parameter 'method'")
}else{
  method <- opt$method
  if(!(method %in% c("median","sum"))){
    stop("The parameter 'method' must be one of 'median','sum'. ")
  }
}

if(file.exists(paste0(outputDir,"/",sampleInfoName))){
  file.remove(paste0(outputDir,"/",sampleInfoName))
  file.copy(from = sampleInfoFile,
            to = paste0(outputDir,"/",sampleInfoName))
}else{
  file.copy(from = sampleInfoFile,
            to = paste0(outputDir,"/",sampleInfoName))
}

data <- readr::read_csv(dataFile,
                        col_types = readr::cols(),
                        progress = FALSE) %>% as.data.frame()
data <- data[!apply(data,1,function(x){all(is.na(x))}),]
data <- data[,!apply(data,2,function(x){all(is.na(x))})]

sampleInfo <-
  readr::read_csv(sampleInfoFile, col_types = readr::cols())%>% as.data.frame()
sampleInfo <- sampleInfo[!apply(sampleInfo,1,function(x){all(is.na(x))}),]
sampleInfo <- sampleInfo[,!apply(sampleInfo,2,function(x){all(is.na(x))})]

# =========================================================================== #
# ============================= data_normalized ============================= #
# =========================================================================== #
#if("QC" %in% unique(sampleInfo$classes)){
#  classes <- c("QC","Sample")
#}else{
#  classes <- c("Sample")
#}
#cal_index <- lapply(classes,
#                    function(x){get_index(colnames(data),sampleInfo,group=x,by="classes")}) %>% unlist
data_info <- data.frame(name=data[,1]) 
if(method == "QC_SVR"){
  # --------------------------- Metnormalizer --------------------------- #
  
  # threads <- parallel::detectCores()[1] - 2
  # Error: segfault from C stack overflow，暂时不能并行，并行会报堆栈溢出的错误。
  # 错误发生的位置/函数： BiocParallel::bplapply
  if(!(all(c("injection.order") %in% colnames(sampleInfo)))){
    stop("Please make sure the column 'injection.order' is in your sample information file.")
  }
  if(!(all(c("compound.id","mz","rt") %in% colnames(data)))){
    stop("Please make sure the column 'compound.id','mz' and 'rt' are all in your data file.")
  }
  threads <- 1
  if("QC" %in% classes){
    # --------------------------- 格式转换 --------------------------- #
    metnor_dir <- paste0(outputDir,"/QC_SVR_input_data")
    if(!dir.exists(metnor_dir)){
      dir.create(metnor_dir,recursive = TRUE)}
    
    new.data <- cbind(name=data$`compound.id`,mz=data$`mz`,rt=data$`rt`,data[,sort(cal_index)])
    sample_index <- get_index(colnames(new.data),sampleInfo,group="Sample",by="classes")
    sample.data <- new.data[,sample_index]
    sample.sd <- apply(sample.data, 1, sd) %>% as.character()
    QC_index <- get_index(colnames(new.data),sampleInfo,group="QC",by="classes")
    QC.data <- new.data[,QC_index]
    QC.sd <- apply(QC.data, 1, sd) %>% as.character()
    sd0 <- union(which(QC.sd==0),which(sample.sd==0))
    if (!sum(sd0==0)) {
       new.data <- new.data
    }else{new.data <- new.data[-sd0,]}
    write.csv(new.data,paste0(metnor_dir,"/data.csv"),row.names = F)
    
    index <- which(sampleInfo$classes %in% classes)
    info <- sampleInfo[index,]
    info$classes[info$classes=="Sample"] <- "Subject"
    info <- info[,colnames(info) %in% c("sample.names","injection.order","classes")]
    colnames(info) <- c("sample.name","injection.order","class")
    #info$sample.name <- make.names(info$sample.name)
    # Note that the data cannot have NAs.
    write.csv(info,paste0(metnor_dir,"/sample.info.csv"),row.names = F)
    
    
    # ---------------------------- metNor ---------------------------- #
    metNor_modified(metnor_dir,threads = threads ,multiple=5,peakplot=peakplot)
    #metNor_modified(metnor_dir,threads = threads ,multiple=1)
    
    # -------------------------- merge data -------------------------- #
    dataFile <- paste0(metnor_dir,"/svr_normalization_result/data_svr_normalization.csv")
    
    data_filted <- readr::read_csv(dataFile,
                                   col_types = readr::cols(),
                                   progress = FALSE) %>% as.data.frame()
    data_filted <- data_filted[,-1]
   
    #if((max(cal_index)+1) < ncol(data)){
    #  data_filted <- cbind(data_filted,data[,(max(cal_index)+1):ncol(data)])
    #}
    
    colnames(data_filted)[c(1,4,5)] <- c("compound.id","sample.rsd","QC.rsd")
  }else{
    stop("QC must be in the colume of 'classes' of sample info file when using method 'QC_SVR' ")
  }
}else{
  switch(method, 
         median = {
           nor_data <- apply(data[,-1],2,
                             FUN = function(x){x/median(x,na.rm = TRUE)})
         },
         sum = {
           nor_data <- apply(data[,-1],2,
                             FUN = function(x){x/sum(x,na.rm = TRUE)})
         })
  
}
data <- cbind(data_info,nor_data)
# ============================================================================ #
# ================================== outputs ================================= #
# ============================================================================ #

write.csv(data,paste0(outputDir,"/data.csv"),row.names = F)


