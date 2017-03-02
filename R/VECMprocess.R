VECMprocess <- 
function(coint.formula, data, stationary.vars = NULL){
  stopifnot(is.ts(data))
  DOLS.list <- DOLSprocess(coint.formula, data)
  DOLS.k <- DOLS.list$DOLS.k
  DOLS.k.HAC <- DOLS.list$DOLS.k.HAC
  y <- DOLS.k$model[,1]
  # construct the yhat from the non-lagged variables and coefficients of DOLS (we ignore nuisance parameters)
  og.Xvars <- which(variable.names(DOLS.k) %in% DOLS.list$data.names)
  yhat <- rowSums(sweep(DOLS.k$model[,og.Xvars], 2, DOLS.k.HAC[og.Xvars,"Estimate"], `*`))
  names(yhat) <- NULL
  yhat.ts <- ts(yhat, start = start(DOLS.k), end = end(DOLS.k), frequency = DOLS.k$frequency)
  # Decompose Error so we can test for asymmetry later
  Error <- y - yhat.ts
  ErrPos <- (diff(Error)>0) * Error
  ErrNeg <- (diff(Error)<=0) * Error 
  Err <- cbind(Error, ErrPos, ErrNeg)
  # lag selection process...
  (maxLags <- floor(dim(DOLS.k$model)[1]^(1/3))) # this will be used for k
  # create the formula dynamically
  if ( !is.null(stationary.vars) ) {
    stationary.vars.vec <- unlist(strsplit(as.character(stationary.vars)[-1]," \\+ "))
  }
  
  ff.maxLags.char <- paste0("diff(",names(DOLS.k$model)[1], ") ~ ", 
                           ifelse("(Intercept)" %in% variable.names(DOLS.k)," 1 + ", "-1 + "),
                           "L(ErrPos,1) + L(ErrNeg,1) + ", 
                           paste0("L(diff(",variable.names(DOLS.k)[og.Xvars],"),1:maxLags)", collapse = " + "),
                           paste0("+ L(diff(",names(DOLS.k$model)[1],"),1:maxLags)"),
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
       VECM.k = VECM.k, # the final model
       VECM.k.HAC = lmtest::coeftest(VECM.k, vcov = sandwich::NeweyWest(VECM.k, lag = k))) # Implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators
}
