#' @title bk_deduct
#' @description Remove background peak using blank.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param data Data matrix.
#' @param blank_index The colume index of blank.
#' @param cal_index The colume index of samples
#' @param fold Cutoff of the sample/blank.
#' @return Clean data after filting.
#' @export
#' @examples
#' bk_deduct(data,blank_index,cal_index,fold=3)

bk_deduct <- function(data,blank_index,cal_index,fold=3){
  before <- nrow(data)
  cat(paste0("Number of original data features: ",before,".\n"))
  data.blank <- data[,blank_index]
  data.cal <- data[,cal_index]
  mean.blank <- apply(data.blank,1,mean,na.rm=TRUE)
  mean.cal <- apply(data.cal,1,mean,na.rm=TRUE)
  FC <- mean.cal/mean.blank
  blank.na.index <- which(apply(data.blank,1,function(x){all(is.na(x))}))
  filted.index <- sort(c(blank.na.index,which(FC > fold)))
  data <- data[filted.index,]
  after <- nrow(data)
  cat(paste0("Number of blank deducted data features: ",after,".\n"))
  cat(paste0(before-after," features have been removed.\n"))
  return(data)
}