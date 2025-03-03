#' @title cor.test_matrix
#' @description get cor and pvalue for two matrix.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param matrixA matrix.
#' @param matrixB matrix.
#' @param alternative one of "two.sided", "less", "greater".
#' @param method one of "pearson", "kendall", "spearman".
#' @return a dataframe contains cor, pvalue, conf.level(only for pearson), method.
#' @export
#' @examples
#' cor.test_matrix(matrixA,alternative = "two.sided",method = "pearson")

cor.test_matrix <- function(matrixA,matrixB=NULL,
                            alternative = c("two.sided", "less", "greater"),
                            method = c("pearson", "kendall", "spearman"))
{
  if(is.null(matrixB)){
    col_names <- colnames(matrixA)
    all_results <- lapply(seq_len(ncol(matrixA)),function(x){
      lapply(seq_len(x), function(y){
        cat(paste0("Calculating ",col_names[x],"_vs_",col_names[y],": method = ",method,", alternative = ",alternative ,"...\n"))
        a <- matrixA[,x]
        b <- matrixA[,y]
        results <- cor.test_modified(a,b,method = method,alternative=alternative)
        results <- as.data.frame(t(results))
        rownames(results) <- paste0(col_names[x],"_vs_",col_names[y])
        return(results)
      })
    })
    all_results <- do.call(rbind,
                           lapply(all_results,function(x){
                             do.call(rbind,x)}))
    return(all_results)
  }else if(!is.null(matrixB)){
    col_namesA <- colnames(matrixA)
    col_namesB <- colnames(matrixB)
    all_results <- lapply(seq_len(ncol(matrixA)),function(x){
      lapply(seq_len(ncol(matrixB)), function(y){
        cat(paste0("Calculating ",col_namesA[x],"_vs_",col_namesB[y],": method = ",method,", alternative = ",alternative ,"...\n"))
        
        a <- matrixA[,x]
        b <- matrixB[,y]
        results <- cor.test_modified(a,b,method = method,alternative=alternative)
        results <- as.data.frame(t(results))
        rownames(results) <- paste0(col_namesA[x],"_vs_",col_namesB[y])
        return(results)
      })
    })
    all_results <- do.call(rbind,
                           lapply(all_results,function(x){
                             do.call(rbind,x)}))
    return(all_results)
  }
}


#' @title cor.test_modified
#' @description get cor and pvalue for two vectors.Base on function of cor.test.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param x vector.
#' @param y vector.
#' @param alternative one of "two.sided", "less", "greater".
#' @param method one of "pearson", "kendall", "spearman".
#' @return a dataframe contains cor, pvalue, conf.level(only for pearson), method.
#' @export
#' @examples
#' cor.test_modified(x,y,alternative = "two.sided",method = "pearson")
cor.test_modified <- function(x,y,
                            alternative = c("two.sided", "less", "greater"),
                            method = c("pearson", "kendall", "spearman"))
{
  if(method=="pearson"){
    results <- cor.test(x,y,method = "pearson",alternative=alternative)
    conf.level <- paste(results$conf.int,collapse = ",")
    results <- c(r=results$estimate,
                          p.values=results$p.value,
                          conf.level=conf.level,
                          method="pearson")
    
  }else if(method=="kendall"){
    results <- cor.test(x,y,method = "kendall",alternative=alternative)
    results <- c(r=results$estimate,
                          p.values=results$p.value,
                          method="kendall")
  }else if(method=="spearman"){
    results <- cor.test(x,y,method = "spearman",alternative=alternative)
    results <- c(r=results$estimate,
                          p.values=results$p.value,
                          method="spearman")
  }
  return(results)
}
