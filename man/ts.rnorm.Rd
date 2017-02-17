\name{ts.rnorm}
\alias{ts.rnorm}
\title{ Time Series Normal Distribution}
\description{
Used the Normal Distribution function to generate random monthly time series datasets that end on the current date.
}
\usage{
 ts.rnorm(n, mean = 0, sd = 1, initial, allownegatives = TRUE)
}
\arguments{
  \item{n}{ number of observations. If length(n) > 1, the length is taken to be the number required. }
  \item{mean}{ vector of means. }
  \item{sd}{ vector of standard deviations. }
  \item{initial}{ numeric; provides the initial value for the time series variable. }
  \item{allownegatives}{  logical; if TRUE (default), negative values remain  otherwise, negative values are replaced with \code{0}. }
}
\details{ If mean or sd are not specified they assume the default values of 0 and 1, respectively. }
\value{ ts.rnorm generates random deviates. }
\seealso{
\code{\link{rnorm} }
}
\examples{
foo1 <- ts.rnorm(n = 150, mean = 0.01, initial = 100)
foo2 <- ts.rnorm(n = 150, mean = 0.01, initial = 283)
foo3 <- ts.rnorm(n = 150, mean = 0.01, initial = 44)
plot(cbind(foo1,foo2,foo3), ylab = "Price", main = "Random Time Series Variables")
}
