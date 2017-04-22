buildDOLS <- 
function (coint.formula, data, robusterrors = FALSE, fixedk = NULL){
  stopifnot(is.ts(data))
  stopifnot(is.null(fixedk)|is.numeric(fixedk))
  output <- list(data.names = colnames(data))
  ff <- coint.formula
  input.vars <- model.frame(formula = ff, data = data)
  RHS <- model.matrix(object = ff, data = input.vars) # just to create the intercept value...if it's there...
  ff.LHS <- colnames(input.vars)[1]
  ff.RHS <- paste(c(ifelse("(Intercept)" %in% colnames(RHS), "1", "-1"), # constant
                    colnames(input.vars[-1]), # input variables
                    paste0("L(diff(", colnames(input.vars)[-1], "),-maxLL:maxLL)")),
                  collapse=" + ")
  ff.maxLL <- paste(ff.LHS, "~", ff.RHS)
  maxLL <- ifelse(is.null(fixedk), floor(dim(input.vars)[1]^(1/3)/2), fixedk) 
  DOLS.maxLL <- dynlm(formula(ff.maxLL), data = data)
  # if user did not enter a fixed number of lags / leads
  if(is.null(fixedk)){
    DOLS.bic <- sapply(1:maxLL, function(maxLL) bic(dynlm(formula(ff.maxLL), data = data, start = start(DOLS.maxLL), end = end(DOLS.maxLL))))
    output$selection <- cbind(`# of lags/leads (k)` = 1:maxLL, SBC = DOLS.bic)
    maxLL.bic <- which.min(DOLS.bic)
    # no need to reestimate, speed optimized here...
    if(maxLL != maxLL.bic){
      output$k <- which.min(maxLL.bic)
      DOLS.k <- dynlm(formula(gsub("maxLL", "maxLL.bic", ff.maxLL)), data = data)
    } else {
      output$k <- maxLL.bic
      DOLS.k <- DOLS.maxLL
    }
  } else {
    output$k <- maxLL
    DOLS.k <- DOLS.maxLL
  }
  output$model <- DOLS.k
  output$model$call <- as.formula(gsub("maxLL", "k", ff.maxLL))
  if(robusterrors){
    output$model$coefficients <- unclass(lmtest::coeftest(DOLS.k, vcov = sandwich::NeweyWest(DOLS.k, lag = output$k))) 
  } else {
    output$model$coefficients <- summary(output$model)$coefficients
  }
  output
}
