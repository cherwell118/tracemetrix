#' @title feature.rm
#' @description Remove features which's NA or 0 beyond the threshold.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param data Data matrix.
#' @param ratio The ratio of NA or 0 in one feature.
#' @param rm.type NA or 0
#' @return data_list.
#' @export
#' @examples
#' feature.rm(data,ratio=0.2,rm.type=0)

feature.rm <- function(data,ratio=0.2,rm.type=0){
  data_list <- list()
  # default: remove variables with missing ratio > 0.2
  if(!(rm.type=="NA"|rm.type==0)){
    stop("rm.type should be NA or 0")
  }else if(rm.type=="NA"){
    ratios <- apply(data,1,function(x){sum(is.na(x))/length(x)})
  }else if(rm.type==0){
    ratios <- apply(data,1,function(x){sum(x==0,na.rm = TRUE)/length(x)})
  }
  cat(paste0("The original feature size is ",nrow(data),".\n"))
  cat(paste0("The original sample size is ",ncol(data),".\n"))
  data_list$original <- data
  data_list$rm.MV <- data[which(!(ratios > ratio)),]
  data_list$row_index <- which(!(ratios > ratio))
  n_rm <- sum(ratios > ratio)
  cat(paste0(n_rm," feature(s) has been eliminated.\n"))
  cat(paste0("The eliminated feature size is ",nrow(data_list$rm.MV),".\n"))
  cat(paste0("The eliminated sample size is ",ncol(data_list$rm.MV),".\n"))
  return(data_list)
}