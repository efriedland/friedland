DOLSprocess <-
function(coint.formula, data) {
  stopifnot(is.ts(data))
  ff <- coint.formula
  input.vars <- model.frame(formula = ff, data = data)
  RHS <- model.matrix(object = ff, data = input.vars) # right hand side
  maxLL <- floor(dim(input.vars)[1]^(1/3)/2) # choosing max lag/lead length
  # maxLL for max Leads/Lags
  # creating the call formula
  ff.maxLL.char <- paste0(colnames(input.vars)[1], " ~ ", # Y variable
                          ifelse("(Intercept)" %in% colnames(RHS)," 1 + ", "-1 + "), # Int
                          paste(colnames(input.vars[-1]), collapse = " + "), # X variables
                          " + ", 
                          paste0("L(diff(", 
                                 colnames(input.vars)[-1],
                                 "),-maxLL:maxLL)", collapse = " + ")
                          ) # Lagged Differences of X variables
  # Saving the DOLS at the max leads & lags so we can fix the start and end date of the selection process
  DOLS.maxLL <- dynlm(as.formula(ff.maxLL.char), data = data)
  # replace "maxLL" with "k"
  ff.k.char <- gsub("maxLL","k",ff.maxLL.char)
  # Retrieve the BIC for each model 
  DOLS.bic <- sapply(1:maxLL, 
                     function(k) bic(dynlm(as.formula(ff.k.char), 
                                           data = data, 
                                           start = start(DOLS.maxLL), 
                                           end = end(DOLS.maxLL))))
  selection <- cbind(`# of lags/leads` = 1:maxLL, SBC = DOLS.bic)
  k <- which.min(DOLS.bic) 
  DOLS.k <- dynlm(formula(ff.k.char), data = data) # best and final model
  list(data.names = colnames(data), # colnames of data to be used for external functions
       call = ff.k.char, # the formula called in the final model
       selection = selection, # BIC results
       leadslags = k, # the selected number of leads/lags
       DOLS.k = DOLS.k, # the final model
       DOLS.k.HAC = lmtest::coeftest(DOLS.k, vcov = sandwich::NeweyWest(DOLS.k, lag = k)) # Implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators
       )
}
