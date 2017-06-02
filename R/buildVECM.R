
buildVECM <- 
  function (coint_formula, data, stationary_vars = NULL, fixedk = NULL, SplitError = TRUE, 
            robusterrors = TRUE, fixedk_DOLS = NULL, selection = AIC) 
{
  stopifnot(is.ts(data))
  stopifnot(is.null(fixedk) || is.numeric(fixedk))
  stopifnot(is.null(fixedk_DOLS) || is.numeric(fixedk_DOLS))
  stopifnot(is.function(selection)) # selection method is a function (should work on a model)
  ff <- coint_formula
  all_names <- dimnames(attr(terms(ff), "factors"))
  y_names <- all_names[[1]][!(all_names[[1]] %in% all_names[[2]])]
  x_names <- all_names[[2]][all_names[[2]] %in% colnames(data)]
  y <- data[, y_names]
  X <- model.matrix(ff, data)
  yhat <- X %*% buildDOLS(coint_formula, data, fixedk = fixedk_DOLS, robusterrors = F)$coefficients[1:dim(X)[2]]
  # Decompose Error term
  Error <- y - ts(yhat, start = start(data), end = end(data), frequency = frequency(data))
  ErrPos <- (diff(y) > 0) * Error
  ErrNeg <- (diff(y) <= 0) * Error
  ff_LHS <- paste0("diff(", y_names, ")")
  ff_RHS <- paste(c(ifelse(attr(terms(ff), "intercept") == 1, "1", "-1"), 
                    ifelse(SplitError, "L(ErrPos, 1) + L(ErrNeg, 1)", "L(Error, 1)"), 
                    paste0("L(diff(", x_names, "), 1:k)", collapse = " + "), 
                    paste0("L(", ff_LHS, ", 1:k)", collapse = " + ")), 
                  collapse = " + ")
  ff_RHS <- ifelse(is.null(stationary_vars), 
                   ff_RHS, 
                   paste0(ff_RHS, " + ", paste0("L(", unlist(strsplit(as.character(stationary_vars)[-1], " \\+ ")), ", 1:k)", 
                                                collapse = " + ")))
  ff_k <- paste(ff_LHS, "~", ff_RHS)
  k <- ifelse(is.null(fixedk), floor(dim(data)[1]^(1/3)), fixedk)
  VECM_k <- dynlm(as.formula(ff_k), data = data)
  VECM_k$Errors <- cbind(Error, ErrPos, ErrNeg)
  if (is.null(fixedk)) {
    k_select <- sapply(1:k, function(k) match.fun(FUN = selection)(dynlm(as.formula(ff_k), data = data, 
                                                                         start = start(VECM_k), end = end(VECM_k))))
    rx <- residuals(VECM_k)
    modselection <- cbind.data.frame(1:k, k_select, VECM_k$df + length(VECM_k$coeff),                                        
                          index2char(index(rx)[1], VECM_k$frequency), 
                          index2char(index(rx)[length(rx)], VECM_k$frequency))
    colnames(modselection) <- c("# of lags (k)", deparse(substitute(selection)), "#Obs", "StartDate", "EndDate")
    if (k != which.min(k_select)) {
      k <- which.min(k_select)
      VECM_k <- dynlm(formula(ff_k), data = data)
    }
  VECM_k$selection <- modselection
  }
  VECM_k$k <- k
  if(robusterrors) VECM_k$HAC <- lmtest::coeftest(VECM_k, vcov = sandwich::NeweyWest(VECM_k, lag = k))
  VECM_k$call <- as.call(c(quote(dynlm), formula = formula(gsub(":k", paste0(":", k), ff_k)), data = substitute(data)))
  class(VECM_k) <- c("workflow", class(VECM_k))
  VECM_k
}
