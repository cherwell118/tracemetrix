#' @title outputFit
#' @description Format the output of model.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param fit Model fitting result
#' @return Data list.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' outputFit(fit)

outputFit  <- function(fit){
  if(class(fit)[1]=="lm"){
    cat("Formating lm")
    # coefficient table
    coefficients <- summary(fit)$coefficients
    conft <- confint(fit)
    if(class(conft) == "matrix"){
        conft <- conft[!apply(conft,1,function(x){all(is.na(x))}),]
    }
    colnames(conft) <- paste(colnames(conft),"Estimate")
    coefficients <- cbind(Estimate=coefficients[,1],conft,coefficients[,-1])
    
    # model related statistics
    ## model p value using F-distribution
    fstatistic <- summary(fit)$fstatistic
    p <- pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=F) %>% signif(.,4)
    
    ## get model expression
    coefs <- coef(fit)[!is.na(coef(fit))] %>% signif(.,3)
    Intercept <- coefs[1]
    coefs <- coefs[-1]
    coefs[coefs > 0] <- paste0("+",coefs[coefs > 0])
    temp <- paste(paste(coefs,names(coefs),sep = "*"),collapse = "")
    equation <- paste0(y," = ",Intercept,temp )
    
    # table of model related statistics
    names <- c("Residual standard error","Multiple R-squared","Adjusted R-squared","F-statistic","p-value","equation")
    values <- c(paste0(signif(summary(fit)$sigma,4)," on ",fit$df.residual," degrees of freedom"),
                signif(summary(fit)$r.squared,4),
                signif(summary(fit)$adj.r.squared,4),
                paste0(signif(fstatistic[1],4)," on ",fstatistic[2]," and ",fstatistic[3], " DF"),
                p,
                equation)
    
    model_stat <- data.frame(names,values)
    
    # return results
    list(coefficients,model_stat)
    
  }else if(class(fit)[1]=="glm"){
    cat("Formating glm")
    coefficients <- summary(fit)$coefficients

    coefs <- coef(fit)[!is.na(coef(fit))]
    OR <- exp(coefs) %>% signif(.,3)

    try_error <- try(confint(fit),silent =TRUE)
    if(!('try-error' %in% class(try_error))){
        conft <- confint(fit)
        if(class(conft) == "matrix"){
          conft <- conft[!apply(conft,1,function(x){all(is.na(x))}),]
        }
        colnames(conft) <- paste(colnames(conft),"Estimate")
        coefficients <- cbind(Estimate=coefficients[,1],conft,coefficients[,-1])

        ORconft <- exp(conft) %>% signif(.,3)
        OR <- cbind(OR,ORconft)
    }
    
    # return results
    results <- list(coefficients,OR)
    names(results) <- c("coefficient","OR")
    return(results)
    
  }else if(class(fit)[1]=="multinom"){
    cat("Formating multinom")
    coefficients <- summary(fit)$coefficients
    Std.Errors <- summary(fit)$standard.errors
    confts <- confint(fit)
    n_levels <- dim(confts)[3]

    results <- lapply((1:n_levels),function(n)multinomFit(n,n_levels,coefficients,Std.Errors,confts))
    
    names(results) <- rownames(coefficients)
    return(results)

  }else if(class(fit)[1]=="polr"){
    cat("Formating polr")
    coefficients <- summary(fit)$coefficients
    p <- pnorm(abs(coefficients[,"t value"]),lower.tail = FALSE)*2
    coefficients <- cbind(coefficients,"Pr(>|t|)"= p)

    coefs <- coef(fit)[!is.na(coef(fit))]
    OR <- exp(coefs) %>% signif(.,3)

    try_error <- try(confint(fit),silent =TRUE)
    if(!('try-error' %in% class(try_error))){
      conft <- confint(fit)
      if(class(conft) == "matrix"){
        conft <- conft[!apply(conft,1,function(x){all(is.na(x))}),]
        colnames(conft) <- paste(colnames(conft),"Estimate")
        coefficients2 <- cbind(Estimate=coefficients[names(coefs),1],conft,coefficients[names(coefs),-1])
      }else{
        names(conft) <- paste(names(conft),"Estimate")
        coefficients2 <- c(Estimate=coefficients[names(coefs),1],conft,coefficients[names(coefs),-1])
      }

    ORconft <- exp(conft) %>% signif(.,3)
    OR <- cbind(OR,ORconft)
}


    # return results    
    if(!('try-error' %in% class(try_error))){
        results <- list(coefficients,coefficients2,OR)
        names(results) <- c("coefficient","coefficients2","OR")
    }else{
        results <- list(coefficients,OR)
        names(results) <- c("coefficient","OR")
    }

    return(results)
  }else{
    stop(paste('Input should be an object of class "lm", "glm","multinom" or "polr".'))
  }
}

#' @title multinomFit
#' @description Format the output of multinom function.
#' @author Ziru Chen
#' \email{chenziru@@picb.ac.cn}
#' @param n Particular level in dependent variant.
#' @param n_levels Number of levels in dependent variant.
#' @param coefficients Results of summary(fit)$coefficients
#' @param Std.Errors Results of summary(fit)$standard.errors
#' @param confts Results of confint(fit)
#' @return Data list.
#' @importFrom magrittr %>%

multinomFit <- function(n,n_levels,coefficients,Std.Errors,confts){
  coefficient <- coefficients[n,]
  name <- rownames(coefficients)[n]
  # conft <- confts[,,n]
  conft <- confts[,,name]
  conft <- conft[!apply(conft,1,function(x){all(is.na(x))}),]
  colnames(conft) <- paste(colnames(conft),"Estimate")
  
  Std.Error <- Std.Errors[name,]
  
  # calculate p
  z <- coefficient/Std.Error
  p <- (1-pnorm(q = abs(z),mean = 0, sd = 1))*2
  
  coefficient_res <- cbind(Estimate=coefficient,conft,"Std. Error"=Std.Error, "z value"=z, "Pr(>|z|)"=p)
  
  OR <- exp(coefficient) %>% signif(.,3)
  ORconft <- exp(conft) %>% signif(.,3)
  OR <- cbind(OR,ORconft)
  
  # return results
  list(coefficient=coefficient_res,OR=OR)
}