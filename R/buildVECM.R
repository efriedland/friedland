buildVECM <- 
  function(coint.formula, data, stationary.vars = NULL, SplitError = FALSE, robusterrors = FALSE, fixedk = NULL){
    stopifnot(is.ts(data))
    stopifnot(is.null(fixedk)|is.numeric(fixedk))
    output <- list(data.names = colnames(data))
    
    DOLS <- buildDOLS(coint.formula, data, robusterrors = F, fixedk)$model
    # save the dependent variables as y
    all.names <- attr(attr(DOLS$terms,"factors"),"dimnames")
    y.names <- all.names[[1]][!(all.names[[1]] %in% all.names[[2]])]
    #y.vars <- which(DOLS.list$data.names %in% y.names)
    y <- data[, y.names] # not y <- DOLS$model[,y.names] as we don't want to start with a limited number of observations from DOLS
    # save the independent variables and multiply them by DOLS coefficients
    og.Xvars <- variable.names(DOLS)[variable.names(DOLS) %in% colnames(data)]
    X <- data[,og.Xvars] 
    
    # to construct the yhat from the non-lagged variables and coefficients of DOLS (we ignore nuisance parameters)
    yhat <- DOLS$coefficients["(Intercept)", "Estimate"] + X %*% DOLS$coefficients[og.Xvars , "Estimate"]
    names(yhat) <- NULL
    yhat.ts <- ts(yhat, start = start(X), end = end(X), frequency = frequency(X))
    # Decompose Error so we can test for asymmetry later
    Error <- y - yhat.ts
    ErrPos <- (diff(y)>0) * Error
    ErrNeg <- (diff(y)<=0) * Error 

    # create the formula dynamically
    if (!is.null(stationary.vars)) {
      stationary.vars.vec <- unlist(strsplit(as.character(stationary.vars)[-1], " \\+ "))
    }
    ff.LHS <- paste0("diff(", names(DOLS$model)[1],")")
    ff.RHS <- paste(c(ifelse("(Intercept)" %in% variable.names(DOLS), "1", "-1"),
                      ifelse(SplitError, "L(ErrorPos, 1) + L(ErrNeg, 1)", "L(Error, 1)"),
                      paste0("L(diff(", og.Xvars, "), 1:maxLags)", collapse = " + "),
                      paste0("L(", ff.LHS, ", 1:maxLags)"),
                      ifelse(is.null(stationary.vars), "", paste0("L(", stationary.vars.vec, ", 1:maxLags)", collapse = ""))), collapse = " + ")
    ff.maxLags <- paste(ff.LHS, "~", ff.RHS)
    maxLags <-  ifelse(is.null(fixedk), floor(dim(DOLS$model)[1]^(1/3)), fixedk)
    VECM.maxLags <- dynlm(as.formula(ff.maxLags), data = data) 
   # if user did not enter a fixed number of lags
    if(is.null(fixedk)){
      VECM.aic <- sapply(1:maxLags, function(maxLags) aic(dynlm(as.formula(ff.maxLags), data = data, start = start(VECM.maxLags), end = end(VECM.maxLags))))
      output$selection <- cbind(`# of lags (k)` = 1:maxLags, AIC = VECM.aic)
      maxLags.aic <- which.min(VECM.aic)
      if(maxLags != maxLags.aic){
        output$k <- maxLags.aic
        VECM.k <- dynlm(formula(gsub("maxLags", "maxLags.aic", ff.maxLags)), data = data)
      } else {
        output$k <- maxLags
        VECM.k <- VECM.maxLags
      }
    } else { # if user did enter a fixed number
      output$k <- maxLags
      VECM.k <- VECM.maxLags
    }
    output$model <- VECM.k
    output$model$call <- as.formula(gsub("maxLags", "k", ff.maxLags))
    if(robusterrors){
      output$model$coefficients <- lmtest::coeftest(VECM.k, vcov = sandwich::NeweyWest(VECM.k, lag = output$k))
    } else {
      output$model$coefficients <- summary(output$model)$coefficients
    }
  output
}
