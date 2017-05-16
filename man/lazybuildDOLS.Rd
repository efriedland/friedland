\name{lazybuildDOLS}
\alias{lazybuildDOLS}
\title{ Comparing DOLS results in the laziest fashion }
\description{ A quick way to compare the DOLS models of each possible cointegrating relationships in your dataset.
The result is a data.frame that is ordered by BIC value.
}
\usage{ lazybuildDOLS(results.coint, dat) }
\arguments{
  \item{results.coint}{ That output of lazyCoint function. }
  \item{dat}{ The time series matrix containing the variables to be used. }
  \item{fixedk}{ The number of lags/leads that the output will have. If is NULL, determines the number of lags through model selection.  }
 }
\details{ 
Relies on buildDOLS and coeftest and NeweyWest HAC estimator to estimate robust errors.
}
\value{
\item{Formula}{ The full formula in character format. }
\item{significance }{ The resulting P-value from the \code{\link{UnitRoot}} function, looks directly to the MacKinnon 1994 table.}
\item{result}{ I(0) == stationary.}
}
\seealso{
\code{\link{buildDOLS}}, \code{\link{lazyCoint}}, \code{\link{bic}}, \code{\link{combn}}
}
\examples{
results <- lazyCoint("logprice", as.ts(M1Germany))
lazybuildDOLS(results, dat = as.ts(M1Germany))
}
