bic <-
function(model){
  n <- df.residual(model) + length(variable.names(model))
  log(sum(resid(model) ^ 2) / n) + length(variable.names(model)) * log(n) / n
}
