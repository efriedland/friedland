buildVECM <- 
function(coint.formula, data, stationary.vars = NULL){
  stopifnot(is.ts(data))
  DOLS.list <- buildDOLS(coint.formula, data)
  DOLS <- DOLS.list$model
  DOLS.HAC <- DOLS.list$robusterrors
  y <- DOLS$model[,1]
  # construct the yhat from the non-lagged variables and coefficients of DOLS (we ignore nuisance parameters)
  og.Xvars <- which(variable.names(DOLS) %in% DOLS.list$data.names)
  if(is.null(dim(DOLS$model[,og.Xvars]))) { # if there is only 1 nonstationary independent variable
    yhat <- DOLS$model[,og.Xvars] * DOLS.HAC[og.Xvars,"Estimate"]
  } else { # if there are more than 1
    yhat <- rowSums(sweep(DOLS$model[,og.Xvars], 2, DOLS.HAC[og.Xvars,"Estimate"], `*`))
  }
  names(yhat) <- NULL
  yhat.ts <- ts(yhat, start = start(DOLS), end = end(DOLS), frequency = DOLS$frequency)
  # Decompose Error so we can test for asymmetry later
  Error <- y - yhat.ts
  ErrPos <- (diff(Error)>0) * Error
  ErrNeg <- (diff(Error)<=0) * Error 
  Err <- cbind(Error, ErrPos, ErrNeg)
  # lag selection process...
  (maxLags <- floor(dim(DOLS$model)[1]^(1/3))) # this will be used for k
  # create the formula dynamically
  if ( !is.null(stationary.vars) ) {
    stationary.vars.vec <- unlist(strsplit(as.character(stationary.vars)[-1]," \\+ "))
  }
  
  ff.maxLags.char <- paste0("diff(",names(DOLS$model)[1], ") ~ ", 
                           ifelse("(Intercept)" %in% variable.names(DOLS)," 1 + ", "-1 + "),
                           "L(ErrPos,1) + L(ErrNeg,1) + ", 
                           paste0("L(diff(",variable.names(DOLS)[og.Xvars],"),1:maxLags)", collapse = " + "),
                           paste0("+ L(diff(",names(DOLS$model)[1],"),1:maxLags)"),
                           ifelse(is.null(stationary.vars),
                                  "",
                                  paste0(" + L(",stationary.vars.vec,",1:maxLags)", collapse = ""))
                           )
  # ifelse() we need the stationary inputs
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
       Error = Err, # A 3 column matrix of the Decomposed Error (y-yhat from the DOLS non-nuissance parameters), only the Positive changes in Error, and the Negative changes in Error
       selection = selection, # AIC results
       lags = k, # the selected number of leads/lags
       model = VECM.k, # the final model
       robusterrors = lmtest::coeftest(VECM.k, vcov = sandwich::NeweyWest(VECM.k, lag = k))) # Implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators
}
