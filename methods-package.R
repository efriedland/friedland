print.workflow <- function(x) {
  if(!is.null(x$selection)){
    cat("\n    Selecting k for Model\n")
    cat("______________________________\n")
    print(x$selection)
  }
  cat("\n    Model with k =",x$k,"\n")
  cat("______________________________")
  print(eval(x$call))
}
                         
summary.workflow <- function(x, ...){
  if("HAC" %in% names(x)){ # if the function has a robusterrors arg
    addHAC <- NextMethod(x) # alter the summary coeffs
    addHAC$coefficients <- x$HAC
    cat("*Summary table depicts HAC estimated errors found by:\n")
    cat(paste0("lmtest::coeftest(model, vcov = sandwich::NeweyWest(model, lag = ",x$k,"))\n"))
    addHAC
  } else {
    warning("\tSummary table does not depict HAC estimated errors\n\tPlease indicate the buildDOLS robusterrors argument to be TRUE")
    print(NextMethod(x))
  }
} 
