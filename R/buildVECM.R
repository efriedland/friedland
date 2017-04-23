buildVECM <- 
  function(coint.formula, data, stationary.vars = NULL, SplitError = FALSE, robusterrors = FALSE, fixedk = NULL){
    stopifnot(is.ts(data))
    stopifnot(is.null(fixedk)|is.numeric(fixedk))
    output <- list(data.names = colnames(data))
    ff <- coint.formula
    # save the dependent variables as y
    all.names <- attr(attr(terms(ff),"factors"),"dimnames")
    y.names <- all.names[[1]][!(all.names[[1]] %in% all.names[[2]])]
    y <- data[, y.names] # not y <- DOLS$model[,y.names] as we don't want to start with a limited number of observations from DOLS
    # save the independent variables and multiply them by DOLS coefficients
    # to construct the yhat from the non-lagged variables and coefficients of DOLS (we ignore nuisance parameters)
    x.names <- all.names[[2]][all.names[[2]] %in% colnames(data)]
    ff <- coint.formula
    X <- model.matrix(ff, data)
    yhat <- X %*% buildDOLS(coint.formula, data, robusterrors = F, fixedk)$model$coefficients[1:dim(X)[2], "Estimate"]
    yhat.ts <- ts(yhat, start = start(data), end = end(data), frequency = frequency(data))
    # Decompose Error so we can test for asymmetry later
    Error <- y - yhat.ts
    ErrPos <- (diff(y)>0) * Error
    ErrNeg <- (diff(y)<=0) * Error 

    # create the formula dynamically
    if (!is.null(stationary.vars)) {
      stationary.vars.vec <- unlist(strsplit(as.character(stationary.vars)[-1], " \\+ "))
    }
    ff.LHS <- paste0("diff(", y.names,")")
    ff.RHS <- paste(c(ifelse(attr(terms(ff), "intercept") == 1, "1", "-1"), 
                  ifelse(SplitError, "L(ErrorPos, 1) + L(ErrNeg, 1)", "L(Error, 1)"), 
                  paste0("L(diff(", x.names, "), 1:maxLags)", collapse = " + "), 
                  paste0("L(", ff.LHS, ", 1:maxLags)", collapse = " + ")), collapse = " + ")
    ff.RHS <- ifelse(is.null(stationary.vars), ff.RHS, paste0(ff.RHS, " + ", paste0("L(", stationary.vars.vec, ", 1:maxLags)", collapse = " + ")))
    ff.maxLags <- paste(ff.LHS, "~", ff.RHS)
    maxLags <-  ifelse(is.null(fixedk), floor(dim(data)[1]^(1/3)), fixedk)
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
