adf <-
function(variable,
                k = 1, 
                drift = TRUE,
                trend = TRUE, 
                startfrom = 0, 
                ADFstatOnly = FALSE){
  stopifnot(is.ts(variable)) # must be time series
  stopifnot(k >= 0) # k must be positive, just k <- abs(k) is not clear enough
  diffvariable <- diff(variable)
  formula <- paste("diffvariable ~ ")
  if(!drift){ formula <- paste(formula, " - 1") }  # If drift = False, then remove the constant
  if(trend){
    s       <- time(diff(variable))
    Time    <- ts(s - s[1], start = s[1], freq = frequency(variable))
    formula <- paste(formula, " + Time")
  }
  formula <- paste(formula, " + L(variable,1)") # 1 Time difference 
  if(k > 0){
    formula <- paste(formula, " + L(diffvariable,1:k)") # Sum of lagged differences
  } 
  x <- dynlm(as.formula(formula), start = startfrom)  
  if(ADFstatOnly){
    ADFstatistic <- summary(x)$coefficients["L(variable, 1)",3]
    ADFstatistic
  } else {
    x
    # remember that the P-Value from traditional summary() is incorrect with adf
  }
}
