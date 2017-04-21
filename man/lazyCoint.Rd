\name{lazyCoint}
\alias{lazyCoint}
\title{ Determining the Cointegrating Relationships in the laziest fashion }
\description{ A quick way to determine all of the possible cointegrating relationships in your dataset, give the dependent variable. 

The function will automatically determine a maximum number of lags based on the number of observations in the dataset. 
It will then use the \code{\link{UnitRootApply}} function (currently set to only include a constant, no trend, using AIC) to determine which variables are nonstationary.
Using the available nonstationary variables, it assembles every possible unique combination of variables to build formulas from.

The result is a matrix that returns only the formulas that appeared to cointegrate with the selected dependent variable. 
}
\usage{ lazyCoint(Y, data) }
\arguments{
  \item{Y}{ A case sensitive character string of the dependent variable. }
  \item{data}{ The time series matrix containing the variables to be used. }
 }
\details{ 
The UnitRootApply and UnitRoot functions used have drift = TRUE and trend = FALSE. 
Relies on MacKinnon 1994 table, which does not specify explicitly what to do when the number of cointegrating variables in the relationship exceeds 6.
}
\value{
\item{Formula}{ The full formula in character format. }
\item{significance }{ The resulting P-value from the \code{\link{UnitRoot}} function, looks directly to the MacKinnon 1994 table.}
\item{result}{ I(0) == stationary.}
}
\seealso{
\code{\link{UnitRootApply}}, \code{\link{UnitRoot}}, \code{\link{MacKinnon}}, \code{\link{combn}}
}
\examples{
lazyCoint("logprice", as.ts(M1Germany))
}
