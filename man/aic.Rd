\name{aic}
\alias{aic}
\alias{bic}
\title{ Information Criterion for model selection }
\description{ Used for model selection. }
\usage{
aic(model)
bic(model)
}
\arguments{
  \item{model}{ A model. }
}
\details{
aic() is to be used when comparing lags. bic() is to be used when comparing lags & leads. 
}
\examples{
## The functions are currently defined as
aic <- function(model){
  n <- df.residual(model) + length(variable.names(model))
  log(sum(resid(model) ^ 2) / n) + (length(variable.names(model)) / n) * 2
} 
bic <- function(model){
    n <- df.residual(model) + length(variable.names(model))
    log(sum(resid(model) ^ 2) / n) + length(variable.names(model)) * log(n) / n
  }
}
