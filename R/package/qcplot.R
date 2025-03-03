#' @title qcplot
#' @description Visualization of data after QC.
#' @author Ziru Chen.
#' \email{chenziru@@picb.ac.cn}
#' @param inputDir Directory of data used for QC plots
#' @param outputDir Output directory for results of QC plots
#' @param dataName Name of data file used for QC plots
#' @param sampleInfoName Name of sample infomation file used for QC plots
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom scales hue_pal
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @examples
#' qcplot(inputDir=outputDir,outputDir=plotDir,dataName="data.csv",sampleInfoName)

qcplot <- function(inputDir,outputDir,dataName,sampleInfoName){
  if(!dir.exists(outputDir)){
    dir.create(outputDir,recursive = TRUE)}
  
  dataFile <- paste0(inputDir,"/",dataName)
  sampleInfoFile <- paste0(inputDir,"/",sampleInfoName)
  
  #=========================================================================#
  #=============================== Read data ===============================#
  #=========================================================================#
  data <- readr::read_csv(dataFile,
                          col_types = readr::cols(),
                          progress = FALSE) %>% as.data.frame()
  data <- data[,!apply(data,2,function(x){all(is.na(x))})]
  data <- data[!apply(data,1,function(x){all(is.na(x))}),]
  
  sampleInfo <-
    readr::read_csv(sampleInfoFile, col_types = readr::cols()) %>%
    dplyr::arrange(injection.order) %>% as.data.frame()
  sampleInfo <- sampleInfo[,!apply(sampleInfo,2,function(x){all(is.na(x))})]
  sampleInfo <- sampleInfo[!apply(sampleInfo,1,function(x){all(is.na(x))}),]
  
  # ======================================================================= #
  # ======================= Setting theme for ggplot2 ===================== #
  # ======================================================================= #
  mytheme <-  theme(legend.position = 'right',
                    #panel.background = element_rect(fill = mycolors[4]),
                    plot.title = element_text(size = rel(2.0), hjust = 0.5),
                    axis.text = element_text(size = rel(1.2)),
                    axis.title = element_text(size = rel(1.4)),
                    legend.title = element_text(size = rel(1.2),hjust = 0),
                    legend.text = element_text(size = rel(1.1)))
  
  #=========================================================================#
  #=============================== RT vs m/z ===============================#
  #=========================================================================#
  rt_mz <- data[,1:3]
  # ---------------------------------------------------------------------- # 
  # use mean intensity of QC as intensity
  # qc_index <- get_index(colnames(data),sampleInfo,group="QC",by="classes")
  # qc <- data[,qc_index]
  # rt_mz$`log10(intensity)` <- apply(qc,1,function(x){mean(x,na.rm = TRUE) %>% log10(.)}) 
  # ---------------------------------------------------------------------- # 
  # use mean intensity of sample as intensity
  sample_index <- get_index(colnames(data),sampleInfo,group="Sample",by="classes")
  sample <- data[,sample_index]
  
  rt_mz$`log10(intensity)` <- apply(sample,1,function(x){mean(x,na.rm = TRUE) %>% log10(.)}) 
  write.csv(rt_mz,paste0(outputDir,"/rt_mz_intensity.csv"),row.names = FALSE)
  
  
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral")[c(9,2)])
  
  pdf(paste0(outputDir,"/rt_mz_intensity.pdf"),width = 12.8,height = 6.2)
  p1 <- ggplot(rt_mz,mapping = aes(x=rt,y = mz)) + 
    geom_point(aes(color = `log10(intensity)`), alpha = 0.7,size = 1.2)
  p2 <- p1 + scale_colour_gradient(low = "grey50", high = "red", na.value = NA)
  p <- p2 + labs(title="Peak profile",y="Mass to charge ratio (m/z)", x = "Retention time") + theme_bw() + mytheme
  print(p)
  dev.off()
  # ========================================================================= #
  # ================================== QC =================================== #
  # ========================================================================= #
  if("QC" %in% unique(sampleInfo$classes)){
    qc_index <- get_index(colnames(data),sampleInfo,group="QC",by="classes")
    qc <- data[,qc_index]
    # -------------------------- intensity boxplot -------------------------- # 
    qc_scale <- scale(t(qc), center = TRUE, 
                      scale = apply(qc, 1, sd, na.rm = TRUE)) %>% t
    
    qc_scale_melt <- reshape2::melt(qc_scale,
                                    value.name = "intensity")
    colnames(qc_scale_melt)[2] <- "QC"
    
    pdf(paste0(outputDir,"/qc_intensity.pdf"))
    p <- ggplot(qc_scale_melt,mapping = aes(x=QC,y = intensity)) + 
      geom_boxplot(color = hue_pal()(4)[1]) + theme_bw()
    p <- p + labs(title="QC intensity boxplot",y="Intensity (auto scaled)", x = "QC") + 
      mytheme +
      theme(axis.text.x = element_text(angle = 45, 
                                       hjust = 0.5, 
                                       vjust = 0.5))
    print(p)
    dev.off()
    
    qc_scale <- data.frame(Compound=data$Compound,qc_scale)
    write.csv(qc_scale,paste0(outputDir,"/qc_intensity_autoScaled.csv"),row.names = FALSE)
    
    # -------------------------- correlation plot --------------------------- # 
    data_cor <- qc
    # data_cor[is.na(data_cor)] <- 0
    qc_cor <-cor(data_cor,method = "pearson",use="complete.obs")
    pdf(paste0(outputDir,"/qc_cor.pdf"))
    pheatmap(qc_cor,cluster_rows = F,cluster_cols = F,
             display_numbers = T,number_format = "%.2f",
             number_color="black")
    dev.off()
    
    write.csv(qc_cor,paste0(outputDir,"/qc_cor.csv"))
  }
  # ======================================================================= #
  # ============================== NA ratios ============================== #
  # ======================================================================= #
  # -------------------------------- peaks -------------------------------- #
  na_ratio_peak <- apply(sample, 1, function(x){(sum(is.na(x))/length(x))*100}) 
  na_ratio_peak <- data.frame(Compound=data$Compound,mv_ratio=na_ratio_peak)
  write.csv(na_ratio_peak,paste0(outputDir,"/missing_value_ratios_peaks.csv"),row.names = FALSE)
  # ------------------------------- samples ------------------------------- # 
  na_ratio_sample <- apply(sample, 2, function(x){(sum(is.na(x))/length(x))*100}) 
  na_ratio_sample <- data.frame(sample = colnames(sample),mv_ratio = na_ratio_sample)
  write.csv(na_ratio_sample,paste0(outputDir,"/missing_value_ratios_sample.csv"),row.names = FALSE)
  # ================================ plots ================================ #
  pdf(paste0(outputDir,"/","missing_value_ratios.pdf"),width = 12.8,height = 6.2)
  # -------------------------------- peaks -------------------------------- #
  p <- ggplot(na_ratio_peak,mapping = aes(x=c(1:nrow(na_ratio_peak)),y = mv_ratio)) + 
    geom_point(color=hue_pal()(4)[1],alpha = 0.6,size = 1.2) + theme_bw()
  p <- p + labs(title="Missing value ratios in peaks",
                y="Missing value ratios in peaks (%)", 
                x = "Peaks") + mytheme
  print(p)
  # ------------------------------- samples ------------------------------- # 
  p <- ggplot(na_ratio_sample,mapping = aes(x=c(1:nrow(na_ratio_sample)),y = mv_ratio)) + 
    geom_point(color=hue_pal()(4)[1],alpha = 0.6,size = 1.2) + theme_bw()
  p <- p + labs(title="Missing value ratios in samples",
                y="Missing value ratios in samples (%)", 
                x = "Samples") + mytheme
  print(p)
  
  dev.off()
  # ======================================================================= #
  # =========================== RSD distribution ========================== #
  # ======================================================================= #
  if("QC" %in% unique(sampleInfo$classes)){
    classes <- c("Sample","QC")
  }else{
    classes <- c("Sample")
  }
  
  all_index <- lapply(classes,
                      function(x){get_index(colnames(data),sampleInfo,group=x,by="classes")}) %>% unlist
  index_data <- data.frame(classes=names(all_index),index=all_index)
  
  data_all <- data[,all_index]
  
  if("QC" %in% unique(sampleInfo$classes)){
    # ---------------------------- samples & QC ----------------------------- # 
    rsd_data <- data.frame(Compound=data[,"Compound"],
                           sample.rsd=data[,"sample.rsd"],
                           QC.rsd=data[,"QC.rsd"])
    # ---------------------------- samples + QC ----------------------------- #
    
    all.rsd <- apply(data_all, 1, rsdFun )
    rsd_data$all.rsd <- all.rsd
    types <- c("sample","QC","all")
  }else{
    # ------------------------------- samples ------------------------------- # 
    rsd_data <- data.frame(Compound=data[,"Compound"],
                           sample.rsd=data[,"sample.rsd"])
    types <- c("sample")
  }
     
    write.csv(rsd_data,paste0(outputDir,"/rsd.csv"),row.names = FALSE)
    
    rsd_cutoff <- 30
    # ------------------------------- rsdPlot ------------------------------- # 
    rsdPlot <- function(type,rsd_cutoff){
      varname <- paste0(type,".pass")
      pass <- ifelse(rsd_data[,paste0(type,".rsd")] < rsd_cutoff,TRUE,FALSE)
      # pass <- ifelse(is.na(rsd_data$rsd),NA,ifelse(rsd_data$rsd < 30,TRUE,FALSE)) %>% as.factor
      passNO <- sum(pass,na.rm = TRUE)
      rsd_data %<>% mutate(., !!varname := pass)
      
      
      p <- ggplot(rsd_data,mapping = aes(x=c(1:nrow(rsd_data)),y = rsd_data[,paste0(type,".rsd")])) + 
        geom_point(aes(color = pass),alpha = 0.5,size = 1.8) + theme_bw() 
      p <- p + geom_hline(aes(yintercept=rsd_cutoff),colour="black", linetype="dashed",size=1.1)
      
      p <- p + labs(title=paste(type,"\n",passNO,"out of",nrow(rsd_data),"peaks with RSDs less than",rsd_cutoff,"%"),
                    y="RSD (%)", 
                    x = "Peaks") + mytheme
      # values=c("#E69F00", "#56B4E9","#999999")
      p <- p + scale_color_manual(values=c(hue_pal()(4)[c(1,3)],"#999999"), 
                                  name="RSD",
                                  breaks=c("FALSE", "TRUE", "NA"),
                                  labels=c(paste0(c("RSD >= ","RSD < "),rsd_cutoff,"%"), "NA"))
      print(p)
    }
    
    pdf(paste0(outputDir,"/RSD.pdf"),width = 12.8,height = 6.2)
    lapply(types, function(type){rsdPlot(type=type,rsd_cutoff=rsd_cutoff)})
    dev.off()
    # ======================================================================= #
    # ================================= PCA ================================= #
    # ======================================================================= #
    # pca_dat <- scale(t(data_all),center=TRUE,scale=TRUE) %>% t
    row.names(data_all) <- data$Compound
    
    pca_dat <- data_all %>% t
    pca_dat <- pca_dat[,!apply(pca_dat,2,function(x){all(is.na(x))})]
    pca_dat <- pca_dat[!apply(pca_dat,1,function(x){all(is.na(x))}),]
    # write.csv(pca_dat,paste0(outputDir,"/data_autoScaled.csv"))
    index_data <- index_data[!apply(pca_dat,1,function(x){all(is.na(x))}),]
    
    naL <- ifelse((is.na(pca_dat) %>% sum) > 1,TRUE,FALSE)
    
    if(naL){
      pca_dat[is.na(pca_dat)] <- 0
      pca_object <- prcomp(pca_dat,center=TRUE,scale=TRUE)
      # 删除缺失值太多的数据
      # naNum <- apply(pca_dat, 2, function(x){(is.na(x) %>% sum)/length(x)})
      # pca_dat <- pca_dat[,(naNum < 0.8)]
      # dim(pca_dat)
      # dim(data_all)
      # test <- nipals(pca_dat)
      # test$c1
      # test <- nipals(pca_dat)
      
      # nb <- estim_ncpPCA(pca_dat,method.cv = "Kfold", verbose = FALSE)
      
      # is.infinite(pca_dat) %>% sum
      # test <- imputePCA(pca_dat,method="EM", ncp=2)
      # test
      # res.pca = FactoMineR::PCA(test$completeObs)
    }else{
      pca_object <- prcomp(pca_dat,center=TRUE,scale=TRUE)
    }
    write.csv(pca_object$x,paste0(outputDir,"/PCA_components_autoScaled.csv"))
    
    pca <- data.frame(pca_object$x[,c(1:2)],class=as.character(index_data$classes))
    write.csv(pca,paste0(outputDir,"/pcaPlotData.csv"))
    
    contr <- summary(pca_object)$importance
    pc1_lab <- paste0("PC1 (",round(contr[2,1]*100,2) ,"%)")
    pc2_lab <- paste0("PC2 (",round(contr[2,2]*100,2) ,"%)")
    
    pc_var <- data.frame(PC=colnames(pca_object$x),Variance=round(contr[2,]*100,2))
    write.csv(pc_var,paste0(outputDir,"/PCs_variance.csv"),row.names = FALSE)
    
    #sample.pcax <- pca[index_data$classes=="Sample",]
    #sample.line <- ellipse::ellipse(cor(sample.pcax$PC1,sample.pcax$PC2),level = 0.95,
    #                                scale = c(sd(sample.pcax$PC1), sd(sample.pcax$PC2)),
    #                                centre = c(mean(sample.pcax$PC1), mean(sample.pcax$PC2))) %>% as.data.frame
    #if("QC" %in% unique(sampleInfo$classes)){
    #  qc.pcax <- pca[index_data$classes=="QC",]
    #  qc.line <- ellipse::ellipse(cor(qc.pcax$PC1,qc.pcax$PC2),level = 0.95,
    #                              scale = c(sd(qc.pcax$PC1), sd(qc.pcax$PC2)),
    #                              centre = c(mean(qc.pcax$PC1), mean(qc.pcax$PC2))) %>% as.data.frame
    #}
    
    
    pdf(paste0(outputDir,"/pca.pdf"))
    p <- ggplot(pca,aes(x=PC1,y=PC2,color = class))+
      theme_bw() +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            strip.background = element_rect(fill = "#0099B47F"),
            strip.text = element_text(color = "white", size = 15))+
      geom_point()+ xlab(pc1_lab)+ylab(pc2_lab)+
      stat_ellipse(show.legend = FALSE,level = 0.95)
    #p + geom_path(data=sample.line,mapping = aes(x = x, y = y), color = 'orange')
    #p + geom_path(data=qc.line,mapping = aes(x = x, y = y), color = 'limegreen')
    # color='dodgerblue'
    print(p)
    dev.off()
  }
  