gatelyPriceAsymm <-
function(x, prefix = "P"){
  stopifnot(is.ts(x)) # must be time series
  Pmax <- 0 # initialize
  for(i in 1:length(x)){ Pmax[i] <- max(x[1:i]) } # records at every date the historical highest price
  Pchange <- diff(Pmax - x) # change in price
  # records price recovery (cumulatively sums every positive change in price)
  Prec    <- cumsum( (Pchange > 0)*Pchange ) 
  # records price cuts (cumulatively sums every negative change in price)
  Pcut    <- cumsum( (Pchange <= 0)*Pchange ) 
  # turns variables into time series objects with the same attributes as the original ts object
  Pmax <- ts(Pmax) ; tsp(Pmax) <- tsp(x) 
  Prec <- ts(Prec) ; tsp(Prec) <- tsp(diff(x))
  Pcut <- ts(Pcut) ; tsp(Pcut) <- tsp(diff(x))
  out <- ts.union(Pmax, Prec, Pcut)
  colnames(out) <- paste0(prefix, c("max","rec","cut"))
  out
}
