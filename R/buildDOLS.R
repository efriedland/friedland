buildDOLS <- 
function (coint.formula, data, robusterrors = FALSE, fixedk = NULL){
  stopifnot(is.ts(data))
  stopifnot(is.null(fixedk)|is.numeric(fixedk))
  output <- list()
  ff <- coint.formula
  all.names <- attr(attr(terms(ff), "factors"), "dimnames")
  y.names <- all.names[[1]][!(all.names[[1]] %in% all.names[[2]])]
  x.names <- all.names[[2]][all.names[[2]] %in% colnames(data)]
  ff.LHS <- y.names
  ff.RHS <- paste(c(ifelse(attr(terms(ff), "intercept") == 1, "1", "-1"), # constant
                    colnames(x.names), # input variables
                    paste0("L(diff(", x.names, "),-maxLL:maxLL)")),
                  collapse=" + ")
  ff.maxLL <- paste(ff.LHS, "~", ff.RHS)
  maxLL <- ifelse(is.null(fixedk), floor(dim(data)[1]^(1/3)/2), fixedk) 
  DOLS.maxLL <- dynlm(formula(ff.maxLL), data = data)
  # if user did not enter a fixed number of lags / leads
  if(is.null(fixedk)){
    DOLS.bic <- sapply(1:maxLL, function(maxLL) bic(dynlm(formula(ff.maxLL), data = data, start = start(DOLS.maxLL), end = end(DOLS.maxLL))))
    output$selection <- cbind(`# of lags/leads (k)` = 1:maxLL, SBC = DOLS.bic)
    maxLL.bic <- which.min(DOLS.bic)
    # no need to reestimate, speed optimized here...
    if(maxLL != maxLL.bic){
      output$k <- maxLL.bic
      DOLS.k <- dynlm(formula(gsub("maxLL", "maxLL.bic", ff.maxLL)), data = data)
    } else { # maxLL == maxLL.bic
      output$k <- maxLL
      DOLS.k <- DOLS.maxLL # no need to re-estimate
    }
  } else { # if user did enter a fixed number
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
