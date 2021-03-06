lazyCoint <-
  function(Y, data, k, ...){
  n <- ifelse(missing(k), floor(nrow(data)^(1/3)), k) # ad hoc method of obtaining max lags to use
  if(missing(k)) warning("No lags supplied (k). Used ",n," lags by default: floor(nrow(",deparse(substitute(data)),")^(1/3))\n")
  Results <- UnitRoot(data, k = n, ...) # determine nonstationary variables
  i1.vars <- names(which(Results[,"result"] == "I(1)"))
  if(!any(i1.vars == Y))
    stop("Your selected Dependent variable (Y) needs to be nonstationary. Please try either ", paste(i1.vars, collapse = " or "))
  X <- i1.vars[-which(i1.vars == Y)]
  possibilities <- unlist(mapply(function(z) {
    combn(X, z, simplify = F) 
  }, 1:length(X)), recursive = F)
  MuMaybe <- matrix(NA, nrow = length(possibilities), ncol = 3)
  Root <- list()
  for(i in 1:length(possibilities)){
    vars <- possibilities[[i]]
    Formula <- formula(paste(Y, "~", paste(vars, collapse = " + ")))
    mu <- residuals(dynlm(Formula, data))
    Root <- UnitRoot(as.ts(mu), k = n, trend = F, cointvariables = min(1 + length(vars),6))
    MuMaybe[i,2] <- round(Root$significance,6)
    MuMaybe[i,3] <- Root$result
    cat(vars,"\n")
  }
  if(is.list(possibilities)){
    MuMaybe[,1] <- paste(Y,"~",rapply(possibilities, paste, collapse = " + "))
  } else {
    MuMaybe[,1] <- paste(Y,"~",possibilities)
  }
  Attempts <- as.data.frame(MuMaybe[which(MuMaybe[,3] == "I(0)"),,drop = F])
  colnames(Attempts) <- c("Formula", "significance", "result")
  Attempts
}
