\name{UnitRoot}
\alias{UnitRoot}
\title{ Unit Root Work Flow}
\description{
UnitRoot requires adf() and MacKinnon and allows the user to quickly determine whether a variable is stationary or non-stationary. The user enters the variable of interest, the model specifications (constant, time trend) and then the number of variables within the relationship. Often 1 variable is used, but ADF tests of residuals will consist of the number of dependent and independent variables in the model. Thus the function can be used to test for cointegration. Calculates ADF statistic for a defined number of models varying by lag number, determines the smallest sized AIC, checks the ADF statistic against MacKinnon Table 3 and determines stationarity. 
}
\usage{
UnitRoot(variable, k, drift = TRUE, trend = TRUE, cointvariables = 1, startfrom = 0)
}
\arguments{
  \item{variable}{ The Time Series variable. }
  \item{k}{ The maximum lags desired, Use only a single positive integer. This will determine the number of models in the selection process. }
  \item{drift}{ Logic statement. If TRUE, includes an intercept/constant in the model. }
  \item{trend}{ Logic statement. If TRUE, includes a time trend in the model. }
  \item{cointvariables}{ Cointegrating variables in the relationship. For the stationarity of a single variable this will be 1. For cointegrating relationships, this equals the number of independent variables + the dependent variable. }
  \item{startfrom}{ The start date of the sample in the test, can be used to limit observations }
}
\details{
Used to pick out which variables are stationary and nonstationary.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{adf}}, \code{\link{aic}}, \code{\link{UnitRootApply}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
