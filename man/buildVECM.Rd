\name{buildVECM}
\alias{buildVECM}
\title{ Building a Vector Error Correction Model }
\description{ A quick way to investigate the short run relationship of cointegrating variables moving to equilibrium. 
             
The function creates the formula format for the VECM equation from the cointegrating relationship. 
Where the cointegrating relationship is \code{Y ~ X1 + X2} and the non-stationary variables are \code{c(z1, z2)},  the VECM equation is \code{diff(Y) ~ 1 + L(ErrPos, 1) + L(ErrNeg, 1) + L(diff(X1), 1:k) + L(diff(X2), 1:k) + L(diff(Y), 1:k) + L(z1, 1:k) + L(z2, 1:k)}, with k = the number of lags.

ErrPos and ErrNeg come from the Error equation Y - Yhat predicted from the coefficients*parameters of a DOLS equation.

It then determines a maximum lag number based of off the number of observations using an ad hoc method and runs the model ad the varying lag sizes, keeping the observation start and end dates restricted to the strictest range. The final number of lags to be chosen is selected by choosing the smallest akaike information criterion value (AIC).

It also calculates the robust standard errors using Newey West with assistance of the lmtest and sandwhich packages. 
}
\usage{ buildVECM(coint_formula, data, stationary_vars = NULL, fixedk = NULL, SplitError = TRUE, robusterrors = TRUE, fixedk_DOLS = NULL, selection = AIC) }
\arguments{
  \item{coint_fomula}{ A representation of the formula that will be built into the VECM short-term equation. }
  \item{data}{ The data frame containing the variables to be called by the formula. }
  \item{stationary_vars}{ A continuation of the formula consisting of stationary variables within data. Must be written with a leading "~" in the format of "~ Z1" to work. }
  \item{fixedk}{The number of lags that the output will have. If is NULL, determines the number of lags through model selection, starting at a maximum value of \code{floor(nrow(data)^(1/3))}.}
  \item{SplitError}{ Default is TRUE. If TRUE will include a positive Error and a negative Error term in the model to highlight asymmetries. }
  \item{robusterrors}{ The HAC estimator.}
  \item{fixedk_DOLS}{The number of lags/leads that the internal DOLS output, which is used to create the error corrrection term, will have. If is NULL, determines the number of lags/leads through model selection, starting at a maximum value of \code{floor(nrow(data)^(1/3)/2)}.}
  \item{selection}{The funciton used during model selection. Default set to BIC but can match any function that acts on a model object}
}
\details{ 
A formula has an implied intercept term. To remove this use either \code{y ~ x - 1} or \code{y ~ 0 + x}.
}
\value{
\item{Error }{returns a 3 column matrix of the decomposed Error (y-yhat from the DOLS non-nuissance parameters), only the Positive changes in Error, and the Negative changes in Error.}
\item{selection }{returns a table showing the lags and the AIC, fixed observations, for choosing the number of lags to use.}
\item{lags }{returns the # of lags identified for the model with the smallest AIC.}
\item{HAC}{returns the robust errors implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators.}
}
\seealso{
\code{\link{buildDOLS}}, \code{\link{aic}}
}
\examples{
dat.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12)
# Use the buildVECM function which does a number of steps displayed in a list object
vecm <- buildVECM(Pgasoline ~ Pcrude, dat.ts, ~ Utilization, SplitError = TRUE, robusterrors = T)
vecm
summary(vecm) 
}
