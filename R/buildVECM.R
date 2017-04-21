buildVECM <- 
  function(coint.formula, data, stationary.vars = NULL, SplitError = FALSE){
    stopifnot(is.ts(data))
    DOLS.list <- buildDOLS(coint.formula, data)
    DOLS <- DOLS.list$model
    DOLS.HAC <- DOLS.list$robusterrors
    # save the dependent variables as y
    all.names <- attr(attr(DOLS$terms,"factors"),"dimnames")
    y.names <- all.names[[1]][!(all.names[[1]] %in% all.names[[2]])]
    #y.vars <- which(DOLS.list$data.names %in% y.names)
    y <- data[, y.names] # not y <- DOLS$model[,y.names] as we don't want to start with a limited number of observations from DOLS
    # save the independent variables and multiply them by DOLS coefficients
    og.Xvars <- variable.names(DOLS)[variable.names(DOLS) %in% colnames(data)]
    X <- data[,og.Xvars] 
    # to construct the yhat from the non-lagged variables and coefficients of DOLS (we ignore nuisance parameters)
    yhat <- DOLS.HAC["(Intercept)", "Estimate"] + X %*% DOLS.HAC[og.Xvars , "Estimate"]
    names(yhat) <- NULL
    yhat.ts <- ts(yhat, start = start(X), end = end(X), frequency = frequency(X))
    # Decompose Error so we can test for asymmetry later
    Error <- y - yhat.ts
    ErrPos <- (diff(y)>0) * Error
    ErrNeg <- (diff(y)<=0) * Error 
    # lag selection process...
    maxLags <- floor(dim(DOLS$model)[1]^(1/3)) # this will be used for k
    # create the formula dynamically
    if ( !is.null(stationary.vars) ) {
      stationary.vars.vec <- unlist(strsplit(as.character(stationary.vars)[-1]," \\+ "))
    }
    ff.maxLags.char <- paste0("diff(",names(DOLS$model)[1], ") ~ ", 
                              ifelse("(Intercept)" %in% variable.names(DOLS)," 1 + ", "-1 + "),
                              ifelse(SplitError, "L(ErrPos,1) + L(ErrNeg,1) + ", "L(Error,1) + "),
                              paste0("L(diff(",og.Xvars,"),1:maxLags)", collapse = " + "),
                              paste0("+ L(diff(",names(DOLS$model)[1],"),1:maxLags)"),
                              ifelse(is.null(stationary.vars),
                                     "",
                                     paste0(" + L(",stationary.vars.vec,",1:maxLags)", collapse = ""))
    )
    # ifelse() lets us alter the formula when it comes to the intercept, the error term, and adding the stationary inputs
    VECM.maxLags <- dynlm(as.formula(ff.maxLags.char), data = data)
    # replace "maxLL" with "k"
    ff.k.char <- gsub("maxLags","k",ff.maxLags.char)
    # Retrieve the BIC for each model 
    VECM.aic <- sapply(1:maxLags, 
                       function(k) aic(dynlm(as.formula(ff.k.char), 
                                             data = data, 
                                             start = start(VECM.maxLags), 
                                             end = end(VECM.maxLags))))
    selection <- cbind(`# of lags` = 1:maxLags, AIC = VECM.aic)
    k <- which.min(VECM.aic) # lags picked
    # run our VECM model with k lags
    VECM.k <- dynlm(formula = as.formula(ff.k.char), data = data)
    
    list(data.names = colnames(data), # colnames of data to be used for external functions
         call = ff.k.char, # the formula called in the final model
         Error = cbind(Error, ErrPos, ErrNeg), # A 3 column matrix of the Decomposed Error (y-yhat from the DOLS non-nuissance parameters), only the Positive changes in Error, and the Negative changes in Error
         selection = selection, # AIC results
         lags = k, # the selected number of leads/lags
         model = VECM.k, # the final model
         robusterrors = lmtest::coeftest(VECM.k, vcov = sandwich::NeweyWest(VECM.k, lag = k))) # Implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators
  }
