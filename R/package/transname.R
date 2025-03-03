#' @title transname 
#' @description Name conversion
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param x data for name conversion
#' @param hashTable table for hash.
#' @return p value.
#' @export
#' @importFrom hash hash
#' @examples
#' transname(x,hashTable)

transname <- function(x,hashTable){
  hashTable <- hash(hashTable[,1],hashTable[,2])
  t <- c()
  for(i in x){
    temp <- hashTable[[i]]
    t <- c(t,temp)
  }
  return(t)
}