\name{buildDOLS}
\alias{buildDOLS}
\title{ Building a Dynamic Ordinary Least Squares Model }
\description{ A quick way to investigate the long run relationship of cointegrating variables using the procedure described by Stock and Watson (1994). Helps determine what the value of the response variables is after it has had enough time to fully adjust to the independent variable. 
             
The function creates the formula format for the DOLS equation from the cointegrating relationship. Where the cointegrating relationship is \code{Y ~ X1 + X2} the DOLS equation is \code{Y ~ 1 + X1 + X2 + L(diff(X1), -k:k) + L(diff(X2), -k:k)}, with k = the number of leads/lags.
             
It then determines a maximum lead/lag number based of off the number of observations using an ad hoc method and runs the model ad the varying lead/lag sizes, keeping the observation start and end dates restricted to the strictest range. The final number of leads/lags to be chosen is selected by choosing the smallest schwarz bayesian information criterion value (SBC/SBIC/SIC).

It also calculates the robust standard errors using Newey West with assistance of the lmtest and sandwhich packages.
} 
\usage{ buildDOLS(coint.formula, data) }
\arguments{
  \item{coint.fomula}{ A symbolic representation of the formula that will be built into the DOLS long-run equation. }
  \item{data}{ The data frame containing the variables to be called by the formula. }
}
\details{ 
A formula has an implied intercept term. To remove this use either \code{y ~ x - 1} or \code{y ~ 0 + x}.
}
\value{
\item{data.names}{returns the column names of the dataset you put in.}
\item{call }{returns the dynlm formula used to make the DOLS.}
\item{selection }{returns a table showing the lags/leads and the SBC, fixed observations, for choosing the number of lags/leads to use.}
\item{leadslags }{returns the # of lags/leads identified for the model with the smallest SBC.}
\item{model }{returns the final model, where you can also find the dates used.}
\item{robusterrors }{returns the robust errors implementing the Newey & West (1987, 1994) heteroskedasticity and autocorrelation consistent (HAC) covariance matrix estimators.}
}
\references{ }
\author{ }
\note{}
\seealso{
\code{\link{buildVECM}}
}
\examples{
data("Asymmetry")
dat.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12)
# Use the DOLS function which does a number of steps
DOLS.list <- buildDOLS(Pgasoline ~ Pcrude, dat.ts)
fit <- DOLS.list$model
start(fit)
end(fit)
summary(fit) # the standard errors are incorrect, instead use
DOLS.list$robusterrors
}
