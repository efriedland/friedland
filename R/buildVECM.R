buildVECM <- 
  function(coint.formula, data, stationary.vars = NULL, fixedk = NULL, SplitError = TRUE, robusterrors = FALSE, fixedk.DOLS = NULL){
    stopifnot(is.ts(data))
    stopifnot(is.null(fixedk)|is.numeric(fixedk))
    ff <- coint.formula
    all.names <- attr(attr(terms(ff),"factors"),"dimnames")
    y.names <- all.names[[1]][!(all.names[[1]] %in% all.names[[2]])]
    x.names <- all.names[[2]][all.names[[2]] %in% colnames(data)]
    y <- data[, y.names] # not y <- DOLS$model[,y.names] as we don't want to start with a limited number of observations from DOLS
    # save the independent variables and multiply them by DOLS coefficients
    # to construct the yhat from the non-lagged variables and coefficients of DOLS (we ignore nuisance parameters)
    X <- model.matrix(ff, data)
    yhat <- X %*% buildDOLS(coint.formula, data, fixedk = fixedk.DOLS, robusterrors = F)$model$coefficients[1:dim(X)[2]]
    yhat.ts <- ts(yhat, start = start(data), end = end(data), frequency = frequency(data))
    # Decompose Error so we can test for asymmetry later
    Error <- y - yhat.ts
    ErrPos <- (diff(y)>0) * Error
    ErrNeg <- (diff(y)<=0) * Error 
    output$Err <- cbind(Error,ErrPos,ErrNeg)
    # create the formula dynamically
    if (!is.null(stationary.vars)) {
      stationary.vars.vec <- unlist(strsplit(as.character(stationary.vars)[-1], " \\+ "))
    }
    ff.LHS <- paste0("diff(", y.names,")")
    ff.RHS <- paste(c(ifelse(attr(terms(ff), "intercept") == 1, "1", "-1"), 
                      ifelse(SplitError, "L(ErrPos, 1) + L(ErrNeg, 1)", "L(Error, 1)"), 
                      paste0("L(diff(", x.names, "), 1:k)", collapse = " + "), 
                      paste0("L(", ff.LHS, ", 1:k)", collapse = " + ")), collapse = " + ")
    ff.RHS <- ifelse(is.null(stationary.vars), ff.RHS, paste0(ff.RHS, " + ", paste0("L(", stationary.vars.vec, ", 1:k)", collapse = " + ")))
    ff.k <- paste(ff.LHS, "~", ff.RHS)
    k <- ifelse(is.null(fixedk), floor(dim(data)[1]^(1/3)), fixedk)
    VECM.k <- dynlm(as.formula(ff.k), data = data) 
    output <- list()
    # Lag Selection method if user did not enter a fixed number
    if(is.null(fixedk)){
      k.aic <- sapply(1:k, function(k) aic(dynlm(as.formula(ff.k), data = data, 
                                                 start = start(VECM.k), 
                                                 end = end(VECM.k))))
      output$selection <- cbind(`# of lags (k)` = 1:k, 
                                AIC = k.aic,
                                `#Obs` = VECM.k$df + length(VECM.k$coeff),
                                StartDate = start(VECM.k)[1],
                                EndDate = end(VECM.k)[1])
      if(k != which.min(k.aic)){
        k <- which.min(k.aic)
        VECM.k <- dynlm(formula(ff.k), data = data)
      }
    }
    output$k <- k
    VECM.k$call <- paste0("dynlm(formula = ",ff.k,", data = ",substitute(data),")")
    output$model <- VECM.k
    if(robusterrors){
      output$robusterrors <- lmtest::coeftest(VECM.k, vcov = sandwich::NeweyWest(VECM.k, lag = output$k))
    }
    output
  }
