suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(mice)))
suppressPackageStartupMessages(suppressWarnings(library(missForest)))
suppressPackageStartupMessages(suppressWarnings(library(imputation)))
suppressPackageStartupMessages(suppressWarnings(library(pcaMethods)))
suppressPackageStartupMessages(suppressWarnings(library(impute)))
suppressPackageStartupMessages(suppressWarnings(library(SeqKnn)))

commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'dataName','d', 1, 'character','Information of data file.',
  'sampleInfoName','a',1,'character','Name of Sample file.',
  'naRatio','r',1,'numeric','NA Ratio cutoff for feature removal',
  'kNum','k',1,'character','parameter in svd,seqknn and knn method.',
  'scaleMethod','s',1,'character','scale method in bpca or ppca,Including "none", "pareto", "vector", "uv"',
  'method','m',1,'character','Method for missing value process. Including "randomForest","mean","median","minimum","halfMinimum","zero"'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

if (!is.null(opt$help)) {
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

if ( is.null(opt$naRatio) ) {
  naRatio <- 0.2
}else{
  naRatio <- opt$naRatio %>% as.numeric
}

if(is.null(opt$method)){
  stop("Please set the parameter 'method'")
}else{
  method <- opt$method
  if(!(method %in% c("randomForest","mean","median","minimum","halfMinimum","zero","KNN","SKNN","SVD","BPCA","PPCA"))){
    stop("The parameter 'method' must be one of randomForest,mean,median,minimum,halfMinimum","zero")
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
  readr::read_csv(sampleInfoFile, col_types = readr::cols()) %>% as.data.frame()
sampleInfo <- sampleInfo[!apply(sampleInfo,1,function(x){all(is.na(x))}),]
sampleInfo <- sampleInfo[,!apply(sampleInfo,2,function(x){all(is.na(x))})]
# ============================================================================ #
# ================================ NA/0 删除 ================================= #
# ============================================================================ #
#if("Sample" %in% unique(sampleInfo$classes)){
#  classes <- "Sample"
#}else{
 # stop("Please check your data. 'Sample' must be contained in colume 'classes'")
#}
#cal_index <- lapply(classes,
#                    function(x){get_index(colnames(data),sampleInfo,group=x,by="classes")}) %>% unlist

data_rm <- feature.rm(data[,-1],ratio=naRatio,rm.type="NA")
data_filted <- data[data_rm$row_index,]

if(is.null(dim(data_filted))){
  stop("Data after filting is empty.")
}

#=====================================================#
#=============== Missing values process ==============#
#=====================================================#
#if(all(c("QC","Sample") %in% unique(sampleInfo$classes))){
#  classes <- c("QC","Sample")
#}else if("Sample" %in% unique(sampleInfo$classes)){
#  classes <- c("Sample")
#}else{
#  stop("Please check your data.")
#}
#cal_index <- lapply(classes,
#                    function(x){get_index(colnames(data_filted),sampleInfo,group=x,by="classes")}) %>% unlist

data_mv <- data_filted[,-1]

# 删除feature(行)全为NA的情况
naRow <- which(apply(data_mv,1,function(x){all(is.na(x))}))
if(length(naRow > 0)){
  data_filted <- data_filted[-c(naRow),]
  data_mv <- data_mv[-c(naRow),]
}
data_info <- data_filted[,1]

Sys.time()
switch(method, 
       #===================== missForest ====================#
       randomForest = {        
         set.seed(123)
         data_mv <- t(data_mv)
         data_mv <- missForest(data_mv,ntree = 100)[[1]] %>% t
       },
       SVD ={
         data_mv <- SVDImpute(data_mv,k=kNum,num.iters = 10,verbose = F)
         data_mv <- data_mv$x
       },       
       PPCA ={
         scaleMethod <- opt$scaleMethod
         if(scaleMethod == "none"){
           data_imp <- pca(data_mv, method="ppca", nPcs=4,scale = scaleMethod,center = FALSE)
           data_mv <- completeObs(data_imp) %>% as.data.frame()
         }else{
           data_imp <- pca(t(data_mv), method= "ppca", nPcs=4,scale = scaleMethod,center = TRUE)
           data_mv <- t(completeObs(data_imp)) %>% as.data.frame()
         }         
       },
       BPCA ={
         scaleMethod <- opt$scaleMethod
         if(scaleMethod == "none"){
           data_imp <- pca(data_mv, method="bpca", nPcs=4,scale = scaleMethod,center = FALSE)
           data_mv <- completeObs(data_imp) %>% as.data.frame()
         }else{
           data_imp <- pca(t(data_mv), method= "bpca", nPcs=4,scale = scaleMethod,center = TRUE)
           data_mv <- t(completeObs(data_imp)) %>% as.data.frame()
         }
       },
#       MI ={
#         names.col <- colnames(data_mv)
#         colnames(data_mv) <- paste0("sample",c(1:ncol(data_mv)))
#         data_imp <- mice(data = data_mv,maxit = 5,m=5,threshold=0.99999999,method = "cart",printFlag = FALSE,seed = 123)
#         data_list <- complete(data_imp,action = "all")
#         data_mv <- data_list[[1]]
#         colnames(data_mv) <- names.col
#         for (i in 1:5) {
#           data_mv_other <- data_list[[i]]
#           colnames(data_mv_other) <- names.col
#           data_mv_other <- cbind(data_info,data_mv_other)
#           write.csv(data_mv_other,paste0(outputDir,"/",data_mv_other,"_MI_datasets_",i,".csv"))
#         }
#       },
       KNN ={
         data_imp <- impute.knn(t(data_mv),k=kNum,colmax = 0.8,rowmax = 0.8,rng.seed =123)
         data_mv <- t(data_imp$data) %>% as.data.frame()
       },
       SKNN ={
         data_mv <- SeqKNN(data_mv,k=kNum)
       },
       mean = {
         data_mv <- apply(data_mv,1,
                          FUN = function(x){ifelse(is.na(x),mean(x,na.rm = TRUE),x)})  %>% t
       },
       median = {
         data_mv <- apply(data_mv,1,
                          FUN = function(x){ifelse(is.na(x),median(x,na.rm = TRUE),x)})  %>% t
       },
       minimum = {
         data_mv <- apply(data_mv,1,
                          FUN = function(x){ifelse(is.na(x),min(x,na.rm = TRUE),x)})  %>% t
       },
       halfMinimum = {
         data_mv <- apply(data_mv,1,
                          FUN = function(x){ifelse(is.na(x),(min(x,na.rm = TRUE) %>% `/`(2)),x)})  %>% t
       },
       zero = {
         data_mv <- apply(data_mv,1,
                          FUN = function(x){ifelse(is.na(x),0,x)}) %>% t
       })
Sys.time()

data <-  cbind(name=data_info,(as.data.frame(data_mv)))
#sampleInfo <- sampleInfo[(sampleInfo$classes %in% classes),]
write.csv(sampleInfo,paste0(outputDir,"/sample_info.csv"),row.names = F)
write.csv(data,paste0(outputDir,"/data.csv"),row.names = F)

                 
                 