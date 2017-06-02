buildDOLS <- 
  function (coint_formula, data, fixedk = NULL, robusterrors = TRUE, selection = BIC){
    stopifnot(is.ts(data)) # time series data
    stopifnot(is.null(fixedk)||is.numeric(fixedk)) # fixed k either is null or is numeric
    stopifnot(is.function(selection)) # selection method is a function (should work on a model)
    # Formula Creation
    ff <- coint_formula 
    all_names <- dimnames(attr(terms(ff), "factors")) # X and Y variables
    y_names <- all_names[[1]][!(all_names[[1]] %in% all_names[[2]])]
    x_names <- all_names[[2]][all_names[[2]] %in% colnames(data)]
    ff_LHS <- y_names
    ff_RHS <- paste(c(ifelse(attr(terms(ff), "intercept") == 1, "1", "-1"), # constant
                      x_names, # input variables
                      paste0("L(diff(", x_names, "),-k:k)")), # lagged differences
                    collapse=" + ")
    ff_k <- paste(ff_LHS, "~", ff_RHS)
    # if k (the maximum number of lags/leads) was not fixed, use a default value
    k <- ifelse(is.null(fixedk), floor(dim(data)[1]^(1/3)/2), fixedk)
    # run the model. If k was fixed, this is the final model:
    DOLS_k <- dynlm(formula(ff_k), data = data)
    # If k was not fixed, DOLS_k will beused to keep constant the start and end dates during model selection
    # Lag/Lead Selection method if used did not enter a fixed number
    if(is.null(fixedk)){
      # Use any selection function that is indicated in the selection argument
      k_select <- sapply(1:k, function(k) match.fun(FUN = selection)(dynlm(formula(ff_k), data = data, 
                                                 start = start(DOLS_k), 
                                                 end = end(DOLS_k))))
      # save the selection matrix results inside the model
      rx <- residuals(DOLS_k)                  
      modselection <- cbind.data.frame(1:k, k_select, DOLS_k$df + length(DOLS_k$coeff), 
                                       index2char(index(rx)[1], DOLS_k$frequency), 
                                       index2char(index(rx)[length(rx)], DOLS_k$frequency))
      colnames(modselection) <- c("# of lags/leads (k)", deparse(substitute(selection)),"#Obs","StartDate", "EndDate")
      # only re-estimate the model if k_select differs from k to be efficient
      if(k != which.min(k_select)){
        k <- which.min(k_select)
        DOLS_k <- dynlm(formula(ff_k), data = data)
      }
      DOLS_k$selection <- modselection                   
    }
  DOLS_k$k <- k # save the lag used inside the model
  # save the HAC estimated errors inside the model
  if(robusterrors) DOLS_k$HAC <- lmtest::coeftest(DOLS_k, vcov = sandwich::NeweyWest(DOLS_k, lag = k))
  # rewriting the call function to be a run on its own 
  DOLS_k$call <- as.call(c(quote(dynlm), formula = formula(gsub("-k:k", paste0("-",k,":",k), ff_k)), data = substitute(data)))                  
  class(DOLS_k) <- c("workflow", class(DOLS_k))
  DOLS_k
}
