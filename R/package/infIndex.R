#' @title infIndex
#' @description Generate data of infIndexPlot.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param fit Name of model.
#' @return Data matrix.
#' @export
#' @examples
#' get_index(fit)

infIndex <- function(fit){
  # cooks.distance
  cooks <- cooks.distance(fit)
  
  # rstudent
  rstu <- rstudent(fit)
  
  # Bonferroni p value
  df <- df.residual(fit) - 1
  rstu <- rstu[!is.na(rstu)]
  n <- length(rstu)
  p <- 2*(pt(abs(rstu), df, lower.tail=FALSE))
  bp <- p.adjust(p,method = "bonferroni")
  
  # hatvalues
  hat_v <- hatvalues(fit)
  
  cbind(cooks.distance=cooks,
             studentized.residuals=rstu,
             p=p,
             "Bonferroni.p-value"=bp,
             "hat-values"=hat_v)
}