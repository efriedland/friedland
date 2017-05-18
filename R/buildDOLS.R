buildDOLS <- 
  function (coint.formula, data, fixedk = NULL, robusterrors = FALSE){
    stopifnot(is.ts(data))
    stopifnot(is.null(fixedk)|is.numeric(fixedk))
    ff <- coint.formula
    all.names <- attr(attr(terms(ff), "factors"), "dimnames")
    y.names <- all.names[[1]][!(all.names[[1]] %in% all.names[[2]])]
    x.names <- all.names[[2]][all.names[[2]] %in% colnames(data)]
    ff.LHS <- y.names
    ff.RHS <- paste(c(ifelse(attr(terms(ff), "intercept") == 1, "1", "-1"), # constant
                      x.names, # input variables
                      paste0("L(diff(", x.names, "),-k:k)")), # lagged differences
                    collapse=" + ")
    ff.k <- paste(ff.LHS, "~", ff.RHS)
    k <- ifelse(is.null(fixedk), floor(dim(data)[1]^(1/3)/2), fixedk)
    DOLS.k <- dynlm(formula(ff.k), data = data) # maxlag or fixedk value 
    output <- list()
    # Lag/Lead Selection method if used did not enter a fixed number
    if(is.null(fixedk)){
      k.bic <- sapply(1:k, function(k) bic(dynlm(formula(ff.k), data = data, 
                                                 start = start(DOLS.k), 
                                                 end = end(DOLS.k))))
      output$selection <- cbind(`# of lags/leads (k)` = 1:k, 
                                SBC = k.bic, 
                                `#Obs` = DOLS.k$df + length(DOLS.k$coeff), 
                                `StartDate` = start(DOLS.k)[1], 
                                `EndDate` = end(DOLS.k)[1])
      if(k != which.min(k.bic)){
        k <- which.min(k.bic)
        DOLS.k <- dynlm(formula(ff.k), data = data)
      } 
    }
    output$k <- k
    DOLS.k$call <- as.call(c(quote(dynlm), 
                             formula = formula(gsub("-k:k", paste0("-",k,":",k), ff.k)), 
                             data = substitute(data)))                  
    output$model <- DOLS.k
    
    if(robusterrors){
      output$robusterrors <- lmtest::coeftest(DOLS.k, vcov = sandwich::NeweyWest(DOLS.k, lag = output$k))
    } 
    output
  }
