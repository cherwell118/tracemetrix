#' @title get_index
#' @description Get index of sample.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param colnames colnames in data file.
#' @param sampleInfo sample information.Include at least two columes: SampleNames, groups
#' @param group Group which you want to get the colume index from data file.
#' @param by Colume of sampleInfo for calulating.
#' @return index numbers.
#' @export
#' @examples
#' get_index(colnames,sampleInfo,group="QC",by="groups")

get_index <- function(colnames,sampleInfo,group="Sample",by="groups"){
  sampleInfo <- as.data.frame(sampleInfo)
  groupInfo <- sampleInfo[sampleInfo[,by]==group,]
  index <- which(colnames %in% groupInfo$sample.names)
  names(index) <- rep(group,length(index))
  return(index)
}