bic <-
function(model){
  n <- df.residual(model) + length(model[["coeff"]])
  log(sum(resid(model) ^ 2) / n) + (length(model[["coeff"]])) * log(n) / n
}
