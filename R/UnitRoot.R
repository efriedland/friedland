UnitRoot <-
function(variable,           # variable
         k,                  # maximum lags
         drift = TRUE,       # constant/intercept
         trend = TRUE,       # time trend
         cointvariables = 1, # number of variables in relationship
         startfrom = 0){
  if(is.null(dim(variable))) {
    n <- ifelse(missing(k), floor(length(variable)^(1/3)), k) # ad hoc method of obtaining max lags to use
    if(missing(k)) warning("No lags supplied (k). Used ", n, "lags by default: floor(length(",deparse(substitute(variable)),")^(1/3))\n")
    k <- n # save myself some time for now.
    M <- matrix(0, ncol=5, nrow = k)
    colnames(M) <- c("Lags", "AIC", "ADF", "#Obs", "Start Date")  
    j <- start(adf(variable, k, drift, trend, startfrom)) 
    for(i in 1:k){
      x <- adf(variable, i, drift, trend, startfrom = j)  
      M[i,"Lags"] <- i
      M[i,"AIC"]  <- aic(x) 
      M[i,"ADF"]  <- adf(variable, i, drift, trend, startfrom = j, ADFstatOnly = TRUE) 
      M[i,"#Obs"] <- x$df + length(x$coeff) # double check why this worked for me
      M[i,"Start Date"] <- time(x)[1]
    }

    Lag <- which.min(M[, "AIC"]) # 
    AIC <- M[Lag, "AIC"] ; names(AIC) <- NULL
    ADFstatistic <- M[which.min( M[, "AIC"] ), "ADF"] ; names(ADFstatistic) <- NULL
    StartDate <- zoo::as.yearmon(M[Lag, "Start Date"]) ; names(StartDate) <- NULL

    # MacKinnon Row lookup type
    if(!drift & !trend) # no constant, no trend
      ltype <- "nc" 
    if(drift) # constant
      ltype <- "c" 
    if(trend) # + trend (will have an error here if the model says drift = F & trend = T)
      ltype <- paste(ltype, "t", sep= "")
    ltype <- paste(ltype, cointvariables,sep = "")
    num <- which(rownames(MacKinnon) == ltype)

    if(ADFstatistic < MacKinnon[num,"Tau_upper"] 
       & ADFstatistic > MacKinnon[num,"Tau_min"]){ # if MacKinnon can be used do the calculation
      h <- MacKinnon[num,1] + MacKinnon[num,2] * ADFstatistic + MacKinnon[num,3] * 0.01 * ADFstatistic ^ 2
      Pvalue = pnorm(h)
      names(Pvalue) <- NULL
      if( pnorm(h) <= 0.05){
        list(selection = M, 
             Lag = Lag,
             StartDate = StartDate,
             AIC = AIC,
             ADFstatistic = ADFstatistic, 
             significance = Pvalue,
             result = "I(0)",
             reason = paste("Using MacKinnon (1994), ",
                            ltype,
                            ", P value of ",
                            format(Pvalue,scientific = TRUE),
                            " is smaller than the 0.05 significance level. Reject Null, variable is I(0)",sep="")) 
      } else {
        list(selection = M, 
             Lag = Lag,
             StartDate = StartDate,
             AIC = AIC,             ADFstatistic = ADFstatistic, 
             significance = Pvalue,
             result = "I(1)",
             reason = paste("Using Mackinnon (1994), ",
                            ltype,
                            ", P value of ",
                            format(Pvalue,scientific = TRUE),
                            " is greater than the 0.05 significance level. Cannot reject Null, variable is I(I)",sep=""))
      }
    } else {
      if(ADFstatistic > MacKinnon[num,"Tau_upper"]){
        list(selection = M, 
             Lag = Lag,
             StartDate = StartDate,
             AIC = AIC,
             ADFstatistic = ADFstatistic, 
             significance = NA,
             result = "I(1)",
             reason = paste("Cannot use MacKinnon (1994), ",
                            ltype,
                            ", ADF Statistic of ",
                            format(ADFstatistic,scientific=TRUE),
                            " too positive, can't reject null, I(1)",sep=""))
      } else{
        list(selection = M, 
             Lag = Lag,
             StartDate = StartDate,
             AIC = AIC,
             ADFstatistic = ADFstatistic, 
             significance = NA,
             result = "I(0)",
             reason = paste("Cannot use MacKinnon (1994), ",
                            ltype,
                            ", ADF Statistic of ",
                            format(ADFstatistic,scientific=TRUE),
                            " too negative, rejects null in favor of alt, I(0)",sep=""))
      }
    }
  } else {
    n <- ifelse(missing(k), floor(nrow(variable)^(1/3)), k) # ad hoc method of obtaining max lags to use
    if(missing(k)) warning("No lags supplied (k). Used ", n, "lags by default: floor(nrow(",deparse(substitute(variable)),")^(1/3))\n")
    t(sapply(lapply(variable, UnitRoot, k = n, drift = drift, trend = trend, 
                    cointvariables = cointvariables, startfrom = startfrom), 
             function(x) x[2:7]))
  }  
}
