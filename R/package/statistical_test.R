#' @title cal_p
#' @description 2 sample statistical test.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param x data of sample 1.
#' @param y data of sample 2.
#' @param paired a logical indicating whether you want a paired t-test
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @return p value.
#' @export
#' @examples
#' cal_p(x,y,paired = FALSE,alternative="two.sided")

cal_p <- function(x,y,paired = FALSE,alternative="two.sided"){
  if(!exists("alternative")){
    #alternative = c("two.sided", "less", "greater")
    alternative="two.sided"
  }
  if(paired){
    shapiro <- shapiro.test((x-y))$p.value
    if(shapiro >= 0.05){
      cat(paste0("t.test(paired=",paired,",alternative = \'",alternative,"\')\n"))
      pvalue <- t.test(x,y,paired = paired,alternative = alternative)$p.value
    }else{
      cat(paste0("wilcox.test(paired=",paired,",alternative = \'",alternative,"\')\n"))
      pvalue <- wilcox.test(x,y,paired = paired,alternative = alternative)$p.value
    }
  }else{
    # 正态性检验
    #Shapiro-Wilk normality test,3<n<5000
    shapiroA <- shapiro.test(x)$p.value
    shapiroB <- shapiro.test(y)$p.value
    if(shapiroA >= 0.05 && shapiroB >= 0.05){
      #F test to compare two variances
      varAB <- var.test(x,y)$p.value
      if(varAB >= 0.05){
        #Two Sample t-test
        cat(paste0("t.test(var.equal = TRUE,paired=",paired,",alternative = \'",alternative,"\')\n"))
        pvalue <- t.test(x,y,var.equal = TRUE,paired = paired,alternative = alternative)$p.value
      }else{
        #Welch Two Sample t-test
        cat(paste0("t.test(var.equal = FALSE,paired=",paired,",alternative = \'",alternative,"\')\n"))
        pvalue <- t.test(x,y,var.equal = FALSE,paired = paired,alternative = alternative)$p.value
      }
    }else{
      cat(paste0("wilcox.test(paired=",paired,",alternative = \'",alternative,"\')\n"))
      pvalue <- wilcox.test(x,y,paired = paired,alternative = alternative)$p.value
    }
  }
  return(pvalue)
}


#' @title get_p
#' @description 2 sample statistical test.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param dataA dataframe 1.
#' @param dataB dataframe 2.
#' @param paired a logical indicating whether you want a paired t-test
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @return p value.
#' @export
#' @examples
#' get_p(dataA,dataB,paired = FALSE,alternative="two.sided")

get_p <- function(dataA,dataB,paired = FALSE,alternative="two.sided"){
  if(nrow(dataA)!=nrow(dataB)){
    stop("The row length of data is not equal.")
  }
  pvalue <- lapply(c(1:nrow(dataA)), function(x){
    cat(paste0(row.names(dataA)[x],": "))
    cal_p(dataA[x,],dataB[x,],paired = paired,alternative=alternative)}) %>% unlist
  return(pvalue)
}

#' @title cal_log2FC
#' @description log2(fold_change) of 2 samples.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param x data of sample 1.
#' @param y data of sample 2.
#' @export
#' @examples
#' cal_log2FC(x,y)

cal_log2FC <- function(x,y){(mean(x)/mean(y)) %>% log2}

#' @title get_log2FC
#' @description log2(fold_change) of 2 dataframes.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param dataA dataframe 1.
#' @param dataB dataframe 2.
#' @return log2FC.
#' @export
#' @examples
#' get_log2FC(x,y)

get_log2FC <- function(dataA,dataB){
  if(nrow(dataA)!=nrow(dataB)){
    stop("The row length of data is not equal.")
  }
  log2FC <- lapply(c(1:nrow(dataA)), function(i){cal_log2FC(dataA[i,],dataB[i,])})%>% unlist
  return(log2FC)
}

#' @title multiTest
#' @description statistics for multiple groups
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param data data.
#' @param groups groups information.
#' @return datalist.
#' @export
#' @examples
#' multiTest(data,groups)

multiTest <- function(data,groups){
  data <- unlist(data)
  if(length(data)!=length(groups)){
    stop("The length of data and groups should be the same.")
  }
  test_data <- data.frame(data=data,groups=groups)
  bart_p <- bartlett.test(data~groups,data = test_data)$p.value
  if(bart_p >= 0.05){
    cat(paste0("use 'aov()' and 'TukeyHSD()'\n"))
    aov_p <- summary(aov(data~groups,data = test_data))[[1]][1,5]
    hsd <- TukeyHSD(aov(data~groups,data = test_data))$groups
    #new_hsd <- hsd %>% t %>% as.vector
    #names(new_hsd) <- paste(rep(row.names(hsd),each=4),colnames(hsd),sep = "_")
    cat(paste0("p = ",aov_p,"\n"))
    cat("Tukey results:\n")
    print(hsd)
    return(list(method="standar anova",p=aov_p,table=hsd))
  }else{
    cat(paste0("use 'oneway.test()' and 'Games-Howell test'\n"))
    oneway_p <- oneway.test(data~groups)$p.value
    GH_test <- tukey(data, groups, method="G")
    Games.Howell <- GH_test$Games.Howell
    group_index <- str_split_fixed(row.names(Games.Howell), ":", 2) %>% as.numeric() # library(stringr)
    gh_name <- levels(groups)[group_index] %>% matrix(.,ncol=2) %>% as.data.frame() %>% dplyr::select("V2","V1")
    gh_name <- do.call(paste,c(gh_name, sep = "-"))
    row.names(Games.Howell) <- gh_name
    GH_test$Games.Howell <- Games.Howell
    cat(paste0("p = ",oneway_p,"\n"))
    cat("Games.Howell results:\n")
    print(Games.Howell)
    cat("\n")
    return(list(method="Welch anova",p=oneway_p,table=GH_test$Games.Howell))
    }
  }