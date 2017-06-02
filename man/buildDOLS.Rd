\name{buildDOLS}
\alias{buildDOLS}
\title{ Building a Dynamic Ordinary Least Squares Model }
\description{ A quick way to investigate the long run relationship of cointegrating variables using the procedure described by Stock and Watson (1994). Helps determine what the value of the response variables is after it has had enough time to fully adjust to the independent variable. 
             
The function creates the formula format for the DOLS equation from the cointegrating relationship. Where the cointegrating relationship is \code{Y ~ X1 + X2} the DOLS equation is \code{Y ~ 1 + X1 + X2 + L(diff(X1), -k:k) + L(diff(X2), -k:k)}, with k = the number of leads/lags.
             
It then determines a maximum lead/lag number based of off the number of observations using an ad hoc method and runs the model ad the varying lead/lag sizes, keeping the observation start and end dates restricted to the strictest range. The final number of leads/lags to be chosen is selected by choosing the smallest schwarz bayesian information criterion value (SBC/SBIC/SIC).

It also calculates the robust standard errors using Newey West with assistance of the lmtest and sandwhich packages.
} 
\usage{ buildDOLS(coint.formula, data, fixedk = NULL, robusterrors = FALSE) }
\arguments{
  \item{coint_fomula}{ A representation of the cointegrating formula that will be built into the DOLS long-run equation. }
  \item{data}{ The data frame containing the variables to be called by the formula. }
  \item{fixedk}{ The number of lags/leads that the output will have. If is NULL, determines the number of lags/leads through model selection, starting at a maximum value of \code{floor(nrow(data)^(1/3)/2)}.. }
  \item{robusterrors}{ The HAC estimator. }
  \item{selection}{ The funciton used during model selection. Default set to BIC but can match any function that acts on a model object.}
}
\details{ 
A formula has an implied intercept term. To remove this use either \code{y ~ x - 1} or \code{y ~ 0 + x}.
}
\value{
\item{selection}{returns a table showing the lags/leads and the selection criteria, fixed observations, and start and end dates used for choosing the number of lags/leads to use.}
\item{k}{returns the # of lags/leads identified for the model with the smallest SBC.}
\item{HAC}{returns the robust errors implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators.}
        }
\seealso{
\code{\link{buildVECM}}, \code{\link{bic}}
}
\examples{
dat.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12)
# Use the DOLS function which does a number of steps
dols <- buildDOLS(Pgasoline ~ Pcrude, dat.ts, robusterrors = TRUE)
dols$selection
summary(dols) # shows HAC estimated errors in the summary output if robusterrors = TRUE
}
