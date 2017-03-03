aic <-
function(model){
  n <- df.residual(model) + length(variable.names(fit$model))
  log(sum(resid(model) ^ 2) / n) + (length(variable.names(fit$model)) / n) * 2
}
