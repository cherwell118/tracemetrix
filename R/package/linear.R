#' @title linear
#' @description Remove background peak using dilution samples.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param x Intensity data of dilution samples.
#' @param digits Digits to keep.
#' @param dilution Dilution concentration gradient.
#' @return Clean data after filting.
#' @export
#' @examples
#' bk_deduct(data,blank_index,cal_index,fold=3)
linear <- function(x,digits=4,dilution=c(0,10,20,50,80,100)){
    if(length(x)==length(dilution)){
        l <- lm(as.numeric(x)~dilution,na.action=na.omit)
        r2 <- summary(l)$r.squared %>% round(.,digits = digits)
        k <- l$coefficients[2] %>%  round(.,digits = digits)
        b <- l$coefficients[1] %>%  round(.,digits = digits)
        result <- cbind(k,b,r2)
        row.names(result) <- ""
        return(result)
    }else{
        stop("The length of x and dilution are inconsistent. Please check your data again.\n")
    }
}